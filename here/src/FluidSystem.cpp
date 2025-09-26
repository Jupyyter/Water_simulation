#include "FluidSystem.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

void advectParticles(std::vector<Particle>& particles, float dt, const std::vector<int>& activeParticles);

FluidSystem::FluidSystem(const Config& config)
    : m_config(config)
{
    m_cellSize = std::max(config.simWidth / 100.0f, config.simHeight / 75.0f);
    m_invCellSize = 1.0f / m_cellSize;
    m_gridWidth = static_cast<int>(std::ceil(config.simWidth * m_invCellSize)) + 1;
    m_gridHeight = static_cast<int>(std::ceil(config.simHeight * m_invCellSize)) + 1;
    int totalCells = m_gridWidth * m_gridHeight;
    m_velocityU.resize(totalCells, 0.0f);
    m_velocityV.resize(totalCells, 0.0f);
    m_prevVelocityU.resize(totalCells, 0.0f);
    m_prevVelocityV.resize(totalCells, 0.0f);
    m_pressure.resize(totalCells, 0.0f);
    m_divergence.resize(totalCells, 0.0f);
    m_particleDensity.resize(totalCells, 0.0f);
    m_cellTypes.resize(totalCells, CellType::AIR);
    m_grid.resize(totalCells, 0.0f);
    m_massU.resize(totalCells, 0.0f);
    m_massV.resize(totalCells, 0.0f);

    m_hashCellSize = config.particleRadius * 2.2f;
    m_hashGridWidth = static_cast<int>(std::ceil(config.simWidth / m_hashCellSize)) + 1;
    m_hashGridHeight = static_cast<int>(std::ceil(config.simHeight / m_hashCellSize)) + 1;
    m_spatialHash.resize(m_hashGridWidth * m_hashGridHeight);
    m_stoneHash.resize(m_hashGridWidth * m_hashGridHeight);

    // Initialize stone grid (each cell can hold one stone particle)
    m_stoneGridCellSize = config.particleRadius * 2.0f; // Exactly particle diameter
    m_stoneGridWidth = static_cast<int>(std::ceil(config.simWidth / m_stoneGridCellSize)) + 1;
    m_stoneGridHeight = static_cast<int>(std::ceil(config.simHeight / m_stoneGridCellSize)) + 1;
    m_stoneGrid.resize(m_stoneGridWidth * m_stoneGridHeight);

    m_particles.reserve(config.maxParticles);
    m_activeParticles.reserve(config.maxParticles);
    m_stoneParticles.reserve(config.maxParticles / 10); // Reserve space for stone particles
    initializeMaterialProperties();
    setupInitialFluid();
    std::cout << "FluidSystem initialized: " << m_gridWidth << "x" << m_gridHeight
        << " grid, cell size: " << m_cellSize 
        << ", stone grid: " << m_stoneGridWidth << "x" << m_stoneGridHeight 
        << ", stone cell size: " << m_stoneGridCellSize << std::endl;
}

FluidSystem::~FluidSystem() {
}

void FluidSystem::initializeMaterialProperties() {
    m_materialProperties.resize(7); // Increased size for STONE
    // Mass values are relative. Sand is 4x heavier than water.
    m_materialProperties[(int)MaterialType::WATER]   = {1.0f, {0.0f, 0.4f, 1.0f}};
    m_materialProperties[(int)MaterialType::SMOKE]   = {0.1f, {0.5f, 0.5f, 0.5f}};
    m_materialProperties[(int)MaterialType::FIRE]    = {0.05f,{1.0f, 0.3f, 0.0f}};
    m_materialProperties[(int)MaterialType::SAND]    = {4.0f, {0.8f, 0.7f, 0.3f}};
    m_materialProperties[(int)MaterialType::OIL]     = {0.8f, {0.2f, 0.1f, 0.0f}};
    m_materialProperties[(int)MaterialType::STONE]   = {10.0f,{0.4f, 0.4f, 0.4f}}; // Heavy, gray stone
    m_materialProperties[(int)MaterialType::TERRAIN] = {10.0f,{0.4f, 0.4f, 0.4f}}; // Same as stone
}

void FluidSystem::setupInitialFluid() {
    // The simulation starts empty.
}

void FluidSystem::update(float deltaTime) {
    if (m_particles.empty()) return;
    deltaTime = std::min(deltaTime, 1.0f / 30.0f);

    m_activeParticles.clear();
    m_stoneParticles.clear();
    
    for (int i = 0; i < static_cast<int>(m_particles.size()); ++i) {
        if (m_particles[i].active) {
            if (m_particles[i].material == MaterialType::STONE || 
                m_particles[i].material == MaterialType::TERRAIN) {
                m_stoneParticles.push_back(i);
            } else {
                m_activeParticles.push_back(i);
            }
        }
    }

    applyGravity(deltaTime);
    advectParticles(m_particles, deltaTime, m_activeParticles); // Only advect non-stone particles
    
    if (m_config.separateParticles) {
        separateParticles(m_config.particleSeparationIterations);
    }

    // Always update stone hash and solid cells to ensure proper collision detection
    buildStoneHash();
    updateSolidCellsFromStones();
    m_stoneHashNeedsUpdate = false; // Reset the flag

    transferVelocitiesToGrid();
    m_prevVelocityU = m_velocityU;
    m_prevVelocityV = m_velocityV;

    enforceGridBoundary();
    solveIncompressibility(m_config.pressureIterations, deltaTime);
    transferVelocitiesFromGrid();

    handleBoundaryCollisions();
    handleStoneCollisions(); // Handle collisions with stone particles
    if (m_config.separateParticles) {
        separateParticles(9);
    }
    updateParticleColors();
    updateSmokeParticles(deltaTime);
    updateFireParticles(deltaTime);
    cleanupInactiveParticles();
}

void FluidSystem::applyGravity(float dt) {
    for (int idx : m_activeParticles) { // Only apply gravity to non-stone particles
        Particle& p = m_particles[idx];
        // Apply gravity scaled by particle mass
        // Heavier particles accelerate faster.
        p.velocity.y += m_config.gravity * p.mass * dt;
        
        // Buoyancy for light materials like smoke/fire
        if (p.material == MaterialType::SMOKE || p.material == MaterialType::FIRE) {
            p.velocity.y += 15.0f * dt;
        }
    }
}

void advectParticles(std::vector<Particle>& particles, float dt, const std::vector<int>& activeParticles) {
    for (int idx : activeParticles) {
        Particle& p = particles[idx];
        p.position += p.velocity * dt;
    }
}

void FluidSystem::separateParticles(int iterations) {
    buildSpatialHash();
    float minDist = 2.0f * m_config.particleRadius;
    float minDist2 = minDist * minDist;

    for (int iter = 0; iter < iterations; ++iter) {
        // Separate active particles from each other
        for (int idx : m_activeParticles) {
            Particle& p1 = m_particles[idx];
            int pxi = static_cast<int>(p1.position.x / m_hashCellSize);
            int pyi = static_cast<int>(p1.position.y / m_hashCellSize);
            
            for (int xi = std::max(pxi - 1, 0); xi <= std::min(pxi + 1, m_hashGridWidth - 1); ++xi) {
                for (int yi = std::max(pyi - 1, 0); yi <= std::min(pyi + 1, m_hashGridHeight - 1); ++yi) {
                    int cellNr = xi * m_hashGridHeight + yi;
                    for (int otherIdx : m_spatialHash[cellNr].particles) {
                        if (idx == otherIdx) continue;
                        Particle& p2 = m_particles[otherIdx];

                        glm::vec2 diff = p1.position - p2.position;
                        float dist2 = glm::dot(diff, diff);

                        if (dist2 > 1e-9 && dist2 < minDist2) {
                            float dist = std::sqrt(dist2);
                            glm::vec2 normal = diff / dist;
                            float displacement_magnitude = 0.5f * (minDist - dist);
                            
                            // Treat stone particles as completely immovable
                            if (p2.material == MaterialType::STONE || p2.material == MaterialType::TERRAIN) {
                                // Stone is immovable - move p1 completely away
                                p1.position += normal * displacement_magnitude * 2.0f;
                                
                                // Extra separation for heavy particles to prevent tunneling
                                if (p1.mass > 1.5f) {
                                    p1.position += normal * m_config.particleRadius * 0.1f;
                                }
                            }
                            else if (p1.material == MaterialType::STONE || p1.material == MaterialType::TERRAIN) {
                                // p1 is stone (immovable) - move p2 completely away
                                p2.position -= normal * displacement_magnitude * 2.0f;
                                
                                // Extra separation for heavy particles
                                if (p2.mass > 1.5f) {
                                    p2.position -= normal * m_config.particleRadius * 0.1f;
                                }
                            }
                            else {
                                // Both are fluid particles - use mass-based separation
                                float total_mass = p1.mass + p2.mass;
                                p1.position += normal * displacement_magnitude * (p2.mass / total_mass);
                                p2.position -= normal * displacement_magnitude * (p1.mass / total_mass);
                            }
                        }
                    }
                }
            }
        }
        
        // Handle collisions between active particles and stone particles
        // using the optimized stone hash for better performance
        for (int activeIdx : m_activeParticles) {
            Particle& activeParticle = m_particles[activeIdx];
            
            // Get the hash cell for this particle
            int pxi = std::clamp(static_cast<int>(activeParticle.position.x / m_hashCellSize), 0, m_hashGridWidth - 1);
            int pyi = std::clamp(static_cast<int>(activeParticle.position.y / m_hashCellSize), 0, m_hashGridHeight - 1);
            
            // Check neighboring cells in stone hash
            for (int xi = std::max(pxi - 1, 0); xi <= std::min(pxi + 1, m_hashGridWidth - 1); ++xi) {
                for (int yi = std::max(pyi - 1, 0); yi <= std::min(pyi + 1, m_hashGridHeight - 1); ++yi) {
                    int cellIdx = xi * m_hashGridHeight + yi;
                    for (int stoneIdx : m_stoneHash[cellIdx].particles) {
                        const Particle& stoneParticle = m_particles[stoneIdx];
                        
                        glm::vec2 diff = activeParticle.position - stoneParticle.position;
                        float dist2 = glm::dot(diff, diff);
                        
                        if (dist2 < minDist2 && dist2 > 1e-9) {
                            float dist = std::sqrt(dist2);
                            glm::vec2 normal = diff / dist;
                            float penetration = minDist - dist;
                            
                            // Move the active particle away from the stone (stone doesn't move)
                            activeParticle.position += normal * penetration;
                            
                            // Add extra separation for heavy particles to prevent tunneling
                            if (activeParticle.mass > 1.5f) {
                                activeParticle.position += normal * m_config.particleRadius * 0.1f;
                            }
                        }
                    }
                }
            }
        }
    }
}

void FluidSystem::transferVelocitiesToGrid() {
    // Reset grids
    std::fill(m_velocityU.begin(), m_velocityU.end(), 0.0f);
    std::fill(m_velocityV.begin(), m_velocityV.end(), 0.0f);
    std::fill(m_massU.begin(), m_massU.end(), 0.0f);
    std::fill(m_massV.begin(), m_massV.end(), 0.0f);

    // Don't reset cell types here - they are now managed by updateSolidCellsFromStones()
    // Only reset FLUID cells to AIR, keep SOLID cells intact
    for (int i = 0; i < static_cast<int>(m_cellTypes.size()); ++i) {
        if (m_cellTypes[i] == CellType::FLUID) {
            m_cellTypes[i] = CellType::AIR;
        }
        // Keep SOLID cells as they are - they represent stone particles
    }

    // Mark fluid cells
    for (int idx : m_activeParticles) {
        const Particle& p = m_particles[idx];
        int xi = std::clamp(static_cast<int>(p.position.x * m_invCellSize), 0, m_gridWidth - 1);
        int yi = std::clamp(static_cast<int>(p.position.y * m_invCellSize), 0, m_gridHeight - 1);
        int cellNr = getGridIndex(xi, yi);
        if (m_cellTypes[cellNr] == CellType::AIR) {
            m_cellTypes[cellNr] = CellType::FLUID;
        }
    }

    // Transfer velocities (only for active/fluid particles)
    for (int idx : m_activeParticles) {
        const Particle& p = m_particles[idx];
        float x = std::clamp(p.position.x, m_cellSize, (m_gridWidth - 1) * m_cellSize);
        float y = std::clamp(p.position.y, m_cellSize, (m_gridHeight - 1) * m_cellSize);

        // Transfer mass-weighted velocity
        transferVelocityComponent(x, y, p.velocity.x, p.mass, m_velocityU, m_massU, 0.0f, m_cellSize * 0.5f);
        transferVelocityComponent(x, y, p.velocity.y, p.mass, m_velocityV, m_massV, m_cellSize * 0.5f, 0.0f);
    }

    for (int i = 0; i < static_cast<int>(m_velocityU.size()); ++i) {
        // Normalize by total mass to get mass-averaged velocity
        if (m_massU[i] > 1e-9) m_velocityU[i] /= m_massU[i];
        if (m_massV[i] > 1e-9) m_velocityV[i] /= m_massV[i];
    }
}
void FluidSystem::updateSolidCellsFromStones() {
    // Only reset SOLID cells to AIR, keep FLUID cells intact
    for (int i = 0; i < static_cast<int>(m_cellTypes.size()); ++i) {
        if (m_cellTypes[i] == CellType::SOLID) {
            m_cellTypes[i] = CellType::AIR;
        }
    }
    
    // Now mark cells as solid where stone particles are located
    for (int idx : m_stoneParticles) {
        const Particle& stone = m_particles[idx];
        int xi = std::clamp(static_cast<int>(stone.position.x * m_invCellSize), 0, m_gridWidth - 1);
        int yi = std::clamp(static_cast<int>(stone.position.y * m_invCellSize), 0, m_gridHeight - 1);
        
        // Mark a small area around the stone particle as solid
        int radius = std::max(1, static_cast<int>(m_config.particleRadius * m_invCellSize));
        for (int dx = -radius; dx <= radius; ++dx) {
            for (int dy = -radius; dy <= radius; ++dy) {
                int x = xi + dx;
                int y = yi + dy;
                if (x >= 0 && x < m_gridWidth && y >= 0 && y < m_gridHeight) {
                    glm::vec2 cellCenter = getGridPosition(x, y) + glm::vec2(m_cellSize * 0.5f);
                    float dist = glm::distance(cellCenter, stone.position);
                    if (dist < m_config.particleRadius * 1.5f) {
                        m_cellTypes[getGridIndex(x, y)] = CellType::SOLID;
                    }
                }
            }
        }
    }
}


void FluidSystem::handleStoneCollisions() {
    // Use optimized collision detection with stone spatial hash
    float restitution = 0.3f;
    float minDist = 2.0f * m_config.particleRadius;
    float minDist2 = minDist * minDist;

    // Check collisions between active particles and stone particles using spatial hash
    for (int activeIdx : m_activeParticles) {
        Particle& activeParticle = m_particles[activeIdx];
        
        // Get the hash cell for this particle
        int pxi = std::clamp(static_cast<int>(activeParticle.position.x / m_hashCellSize), 0, m_hashGridWidth - 1);
        int pyi = std::clamp(static_cast<int>(activeParticle.position.y / m_hashCellSize), 0, m_hashGridHeight - 1);
        
        // Check neighboring cells in stone hash
        for (int xi = std::max(pxi - 1, 0); xi <= std::min(pxi + 1, m_hashGridWidth - 1); ++xi) {
            for (int yi = std::max(pyi - 1, 0); yi <= std::min(pyi + 1, m_hashGridHeight - 1); ++yi) {
                int cellIdx = xi * m_hashGridHeight + yi;
                for (int stoneIdx : m_stoneHash[cellIdx].particles) {
                    const Particle& stoneParticle = m_particles[stoneIdx];
                    
                    glm::vec2 diff = activeParticle.position - stoneParticle.position;
                    float dist2 = glm::dot(diff, diff);
                    
                    if (dist2 < minDist2 && dist2 > 1e-9) {
                        float dist = std::sqrt(dist2);
                        glm::vec2 normal = diff / dist;
                        float penetration = minDist - dist;
                        
                        // Move the active particle away from the stone
                        activeParticle.position += normal * penetration;
                        
                        // Apply collision response (stone doesn't move)
                        float vn = glm::dot(activeParticle.velocity, normal);
                        if (vn < 0) {
                            activeParticle.velocity -= (1.0f + restitution) * vn * normal;
                        }
                    }
                }
            }
        }
    }
}

void FluidSystem::transferVelocityComponent(float x, float y, float vel, float mass,
    std::vector<float>& field, std::vector<float>& mass_field,
    float offsetX, float offsetY) {
    float fx = (x - offsetX) * m_invCellSize;
    float fy = (y - offsetY) * m_invCellSize;

    int x0 = std::min(static_cast<int>(fx), m_gridWidth - 2);
    int y0 = std::min(static_cast<int>(fy), m_gridHeight - 2);
    float tx = fx - x0; float ty = fy - y0;
    float sx = 1.0f - tx; float sy = 1.0f - ty;
    float d0 = sx * sy; float d1 = tx * sy;
    float d2 = tx * ty; float d3 = sx * ty;

    int nr0 = getGridIndex(x0, y0);
    int nr1 = getGridIndex(x0 + 1, y0);
    int nr2 = getGridIndex(x0 + 1, y0 + 1);
    int nr3 = getGridIndex(x0, y0 + 1);

    // Add mass-weighted velocity and mass to the grid
    field[nr0] += vel * mass * d0; mass_field[nr0] += mass * d0;
    field[nr1] += vel * mass * d1; mass_field[nr1] += mass * d1;
    field[nr2] += vel * mass * d2; mass_field[nr2] += mass * d2;
    field[nr3] += vel * mass * d3; mass_field[nr3] += mass * d3;
}

void FluidSystem::transferVelocitiesFromGrid() {
    // Only transfer velocities back to active (non-stone) particles
    for (int idx : m_activeParticles) {
        Particle& p = m_particles[idx];

        float x = std::clamp(p.position.x, m_cellSize, (m_gridWidth - 2) * m_cellSize);
        float y = std::clamp(p.position.y, m_cellSize, (m_gridHeight - 2) * m_cellSize);

        float picU = interpolateVelocityComponent(x, y, m_velocityU, 0.0f, m_cellSize * 0.5f);
        float picV = interpolateVelocityComponent(x, y, m_velocityV, m_cellSize * 0.5f, 0.0f);
        float prevU = interpolateVelocityComponent(x, y, m_prevVelocityU, 0.0f, m_cellSize * 0.5f);
        float prevV = interpolateVelocityComponent(x, y, m_prevVelocityV, m_cellSize * 0.5f, 0.0f);

        float deltaU = picU - prevU;
        float deltaV = picV - prevV;
        float flipU = p.velocity.x + deltaU;
        float flipV = p.velocity.y + deltaV;

        p.velocity.x = (1.0f - m_config.flipRatio) * picU + m_config.flipRatio * flipU;
        p.velocity.y = (1.0f - m_config.flipRatio) * picV + m_config.flipRatio * flipV;
    }
}

void FluidSystem::addParticles(const glm::vec2& position, float radius, int count, MaterialType material) {
    // Treat TERRAIN as STONE
    if (material == MaterialType::TERRAIN) {
        material = MaterialType::STONE;
    }

    // Special handling for stone particles - use grid-based placement
    if (material == MaterialType::STONE) {
        // Find all empty stone grid cells within the spawn radius
        glm::ivec2 centerGrid = worldToStoneGrid(position);
        int gridRadius = static_cast<int>(std::ceil(radius / m_stoneGridCellSize));
        
        std::vector<glm::ivec2> emptyCells;
        
        // Collect all empty cells within the circular area
        for (int dx = -gridRadius; dx <= gridRadius; ++dx) {
            for (int dy = -gridRadius; dy <= gridRadius; ++dy) {
                glm::ivec2 gridPos = centerGrid + glm::ivec2(dx, dy);
                
                // Check if within bounds
                if (gridPos.x >= 0 && gridPos.x < m_stoneGridWidth && 
                    gridPos.y >= 0 && gridPos.y < m_stoneGridHeight) {
                    
                    // Check if within circular radius
                    glm::vec2 cellWorldPos = stoneGridToWorld(gridPos);
                    float dist = glm::distance(cellWorldPos, position);
                    
                    if (dist <= radius && !isStoneGridCellOccupied(gridPos.x, gridPos.y)) {
                        emptyCells.push_back(gridPos);
                    }
                }
            }
        }
        
        // Create stone particles in empty cells (up to the requested count)
        int stonesToCreate = std::min(count, static_cast<int>(emptyCells.size()));
        stonesToCreate = std::min(stonesToCreate, static_cast<int>(m_config.maxParticles - m_particles.size()));
        
        for (int i = 0; i < stonesToCreate; ++i) {
            glm::ivec2 gridPos = emptyCells[i];
            glm::vec2 stonePos = stoneGridToWorld(gridPos);
            
            Particle p;
            p.position = stonePos;
            p.velocity = glm::vec2(0.0f);
            p.material = material;
            p.mass = m_materialProperties[static_cast<int>(material)].mass;
            p.color = m_materialProperties[static_cast<int>(material)].color;
            p.active = true;
            p.life = 1.0f;
            p.isStatic = true;

            int particleIndex = m_particles.size();
            m_particles.push_back(p);
            
            // Mark the grid cell as occupied
            addStoneToGrid(gridPos.x, gridPos.y, particleIndex);
        }
        
        // Flag that stone hash needs updating
        if (stonesToCreate > 0) {
            m_stoneHashNeedsUpdate = true;
        }
        
        return;
    }
    
    // Original particle spawning for non-stone materials
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-1.0f, 1.0f);

    for (int i = 0; i < count && m_particles.size() < m_config.maxParticles; ++i) {
        Particle p;

        float angle = dis(gen) * 3.14159f * 2.0f;
        float r = std::sqrt(std::abs(dis(gen))) * radius;
        p.position = position + glm::vec2(std::cos(angle) * r, std::sin(angle) * r);
        p.position.x = std::clamp(p.position.x, m_config.particleRadius, m_config.simWidth - m_config.particleRadius);
        p.position.y = std::clamp(p.position.y, m_config.particleRadius, m_config.simHeight - m_config.particleRadius);

        p.velocity = glm::vec2(0.0f);
        p.material = material;
        p.mass = m_materialProperties[static_cast<int>(material)].mass;
        p.color = m_materialProperties[static_cast<int>(material)].color;
        p.active = true;
        p.life = 1.0f;
        p.isStatic = false;

        m_particles.push_back(p);
    }
}

// Backwards compatibility functions
void FluidSystem::addStaticTerrain(const glm::vec2& position, float radius) {
    // Convert to stone particles
    int stoneCount = static_cast<int>(radius / (m_config.particleRadius * 1.5f)) * 8; // More stones for larger radius
    addParticles(position, radius, stoneCount, MaterialType::STONE);
}

void FluidSystem::removeTerrainInArea(const glm::vec2& position, float radius) {
    // Remove stone particles in area
    removeParticlesInArea(position, radius);
}

// Rest of the functions remain mostly the same...

void FluidSystem::cleanupInactiveParticles() {
    m_particles.erase(
        std::remove_if(m_particles.begin(), m_particles.end(),
            [](const Particle& p) { return !p.active; }),
        m_particles.end()
    );
}

void FluidSystem::enforceGridBoundary() {
    for (int i = 0; i < m_gridWidth; ++i) {
        for (int j = 0; j < m_gridHeight; ++j) {
            if (m_cellTypes[getGridIndex(i,j)] == CellType::SOLID) {
                m_velocityU[getGridIndex(i,j)] = 0.0f;
                m_velocityV[getGridIndex(i,j)] = 0.0f;
                if (i + 1 < m_gridWidth) m_velocityU[getGridIndex(i + 1, j)] = 0.0f;
                if (j + 1 < m_gridHeight) m_velocityV[getGridIndex(i, j + 1)] = 0.0f;
            }
        }
    }
}

void FluidSystem::solveIncompressibility(int iterations, float dt) {
    std::fill(m_pressure.begin(), m_pressure.end(), 0.0f);
    for (int iter = 0; iter < iterations; ++iter) {
        for (int i = 1; i < m_gridWidth - 1; ++i) {
            for (int j = 1; j < m_gridHeight - 1; ++j) {
                if (m_cellTypes[getGridIndex(i,j)] != CellType::FLUID) continue;
                int center = getGridIndex(i, j);
                int left = getGridIndex(i - 1, j);
                int right = getGridIndex(i + 1, j);
                int bottom = getGridIndex(i, j - 1);
                int top = getGridIndex(i, j + 1);
                float sx0 = (m_cellTypes[left] != CellType::SOLID);
                float sx1 = (m_cellTypes[right] != CellType::SOLID);
                float sy0 = (m_cellTypes[bottom] != CellType::SOLID);
                float sy1 = (m_cellTypes[top] != CellType::SOLID);
                float s = sx0 + sx1 + sy0 + sy1;
                if (s == 0.0f) continue;
                float divergence = m_velocityU[right] - m_velocityU[center] + m_velocityV[top] - m_velocityV[center];
                float p = -divergence / s;
                p *= m_overRelaxation;
                m_velocityU[center] -= sx0 * p;
                m_velocityU[right] += sx1 * p;
                m_velocityV[center] -= sy0 * p;
                m_velocityV[top] += sy1 * p;
            }
        }
    }
}

void FluidSystem::handleBoundaryCollisions() {
    float margin = m_config.particleRadius;
    float restitution = -m_config.particleRestitution;
    for (int idx : m_activeParticles) {
        Particle& p = m_particles[idx];
        if (p.position.x < margin) {
            p.position.x = margin;
            if (p.velocity.x < 0) p.velocity.x *= restitution;
        } else if (p.position.x > m_config.simWidth - margin) {
            p.position.x = m_config.simWidth - margin;
            if (p.velocity.x > 0) p.velocity.x *= restitution;
        }
        if (p.position.y < margin) {
            p.position.y = margin;
            if (p.velocity.y < 0) p.velocity.y *= restitution;
        } else if (p.position.y > m_config.simHeight - margin) {
            p.position.y = m_config.simHeight - margin;
            if (p.velocity.y > 0) p.velocity.y *= restitution;
        }
    }
}

void FluidSystem::buildSpatialHash() {
    for (auto& cell : m_spatialHash) {
        cell.particles.clear();
    }
    // Only add fluid particles to the spatial hash
    for (int idx : m_activeParticles) {
        const Particle& p = m_particles[idx];
        int hashX = std::clamp(static_cast<int>(p.position.x / m_hashCellSize), 0, m_hashGridWidth - 1);
        int hashY = std::clamp(static_cast<int>(p.position.y / m_hashCellSize), 0, m_hashGridHeight - 1);
        int hashIndex = hashX * m_hashGridHeight + hashY;
        m_spatialHash[hashIndex].particles.push_back(idx);
    }
}

void FluidSystem::buildStoneHash() {
    for (auto& cell : m_stoneHash) {
        cell.particles.clear();
    }
    // Only add stone particles to the stone hash
    for (int idx : m_stoneParticles) {
        const Particle& p = m_particles[idx];
        int hashX = std::clamp(static_cast<int>(p.position.x / m_hashCellSize), 0, m_hashGridWidth - 1);
        int hashY = std::clamp(static_cast<int>(p.position.y / m_hashCellSize), 0, m_hashGridHeight - 1);
        int hashIndex = hashX * m_hashGridHeight + hashY;
        m_stoneHash[hashIndex].particles.push_back(idx);
    }
}

void FluidSystem::updateParticleColors() {
    for (int i = 0; i < static_cast<int>(m_particles.size()); ++i) {
        if (!m_particles[i].active) continue;
        Particle& p = m_particles[i];
        if (p.material != MaterialType::FIRE) {
            p.color = m_materialProperties[static_cast<int>(p.material)].color;
        }
    }
}

void FluidSystem::updateSmokeParticles(float dt) {
    for (int idx : m_activeParticles) {
        Particle& p = m_particles[idx];
        if (p.material == MaterialType::SMOKE) {
            p.life -= 0.5f * dt;
            if (p.life < 0.0f) p.active = false;
        }
    }
}

void FluidSystem::updateFireParticles(float dt) {
    for (int idx : m_activeParticles) {
        Particle& p = m_particles[idx];
        if (p.material == MaterialType::FIRE) {
            p.life -= (1.5f * (1.1f - p.life)) * dt;
            if (p.life <= 0.0f) {
                p.material = MaterialType::SMOKE;
                p.life = 1.0f;
                p.color = m_materialProperties[static_cast<int>(MaterialType::SMOKE)].color;
            } else {
                float intensity = std::max(0.0f, p.life);
                p.color = glm::vec3(1.0f, 0.3f + 0.4f * intensity, 0.0f);
            }
        }
    }
}

void FluidSystem::removeParticlesInArea(const glm::vec2& position, float radius) {
    float radius2 = radius * radius;
    bool removedStone = false;
    
    for (int i = 0; i < static_cast<int>(m_particles.size()); ++i) {
        Particle& particle = m_particles[i];
        if (particle.active) {
            glm::vec2 diff = particle.position - position;
            if (glm::dot(diff, diff) < radius2) {
                if (particle.material == MaterialType::STONE || particle.material == MaterialType::TERRAIN) {
                    // Remove from stone grid
                    glm::ivec2 gridPos = worldToStoneGrid(particle.position);
                    if (gridPos.x >= 0 && gridPos.x < m_stoneGridWidth && 
                        gridPos.y >= 0 && gridPos.y < m_stoneGridHeight) {
                        removeStoneFromGrid(gridPos.x, gridPos.y);
                    }
                    removedStone = true;
                }
                particle.active = false;
            }
        }
    }
}

int FluidSystem::getParticleCount() const {
    return m_particles.size();
}

float FluidSystem::interpolateVelocityComponent(float x, float y, const std::vector<float>& field, float offsetX, float offsetY) const {
    float fx = (x - offsetX) * m_invCellSize;
    float fy = (y - offsetY) * m_invCellSize;
    int x0 = std::min(static_cast<int>(fx), m_gridWidth - 2);
    int y0 = std::min(static_cast<int>(fy), m_gridHeight - 2);
    int x1 = std::min(x0 + 1, m_gridWidth - 2);
    int y1 = std::min(y0 + 1, m_gridHeight - 2);
    float tx = fx - x0;
    float ty = fy - y0;
    float sx = 1.0f - tx;
    float sy = 1.0f - ty;
    float v00 = field[getGridIndex(x0, y0)];
    float v10 = field[getGridIndex(x1, y0)];
    float v01 = field[getGridIndex(x0, y1)];
    float v11 = field[getGridIndex(x1, y1)];
    return sx * sy * v00 + tx * sy * v10 + tx * ty * v11 + sx * ty * v01;
}

glm::vec2 FluidSystem::getGridPosition(int x, int y) const {
    return glm::vec2(x * m_cellSize, y * m_cellSize);
}

glm::ivec2 FluidSystem::worldToGrid(const glm::vec2& worldPos) const {
    return glm::ivec2(static_cast<int>(worldPos.x * m_invCellSize), static_cast<int>(worldPos.y * m_invCellSize));
}

glm::vec2 FluidSystem::gridToWorld(const glm::ivec2& gridPos) const {
    return glm::vec2(gridPos) * m_cellSize;
}

void FluidSystem::reset() {
    // Clear all particles
    m_particles.clear();
    m_activeParticles.clear();
    m_stoneParticles.clear();
    
    // Clear stone grid
    for (auto& cell : m_stoneGrid) {
        cell.occupied = false;
        cell.particleIndex = -1;
    }
    
    // Flag that stone hash needs updating
    m_stoneHashNeedsUpdate = true;

    // Reset all simulation grids to their initial state
    std::fill(m_velocityU.begin(), m_velocityU.end(), 0.0f);
    std::fill(m_velocityV.begin(), m_velocityV.end(), 0.0f);
    std::fill(m_pressure.begin(), m_pressure.end(), 0.0f);
    std::fill(m_particleDensity.begin(), m_particleDensity.end(), 0.0f);
    std::fill(m_cellTypes.begin(), m_cellTypes.end(), CellType::AIR);

    // Reset simulation state variables
    m_particleRestDensity = 0.0f;

    // Call the setup function (which is currently empty)
    setupInitialFluid();
}

// Stone grid helper functions
glm::ivec2 FluidSystem::worldToStoneGrid(const glm::vec2& worldPos) const {
    return glm::ivec2(
        static_cast<int>(worldPos.x / m_stoneGridCellSize),
        static_cast<int>(worldPos.y / m_stoneGridCellSize)
    );
}

glm::vec2 FluidSystem::stoneGridToWorld(const glm::ivec2& gridPos) const {
    return glm::vec2(
        (gridPos.x + 0.5f) * m_stoneGridCellSize,
        (gridPos.y + 0.5f) * m_stoneGridCellSize
    );
}

bool FluidSystem::isStoneGridCellOccupied(int x, int y) const {
    if (x < 0 || x >= m_stoneGridWidth || y < 0 || y >= m_stoneGridHeight) {
        return true; // Consider out-of-bounds as occupied
    }
    return m_stoneGrid[getStoneGridIndex(x, y)].occupied;
}

void FluidSystem::addStoneToGrid(int x, int y, int particleIndex) {
    if (x >= 0 && x < m_stoneGridWidth && y >= 0 && y < m_stoneGridHeight) {
        int index = getStoneGridIndex(x, y);
        m_stoneGrid[index].occupied = true;
        m_stoneGrid[index].particleIndex = particleIndex;
    }
}

void FluidSystem::removeStoneFromGrid(int x, int y) {
    if (x >= 0 && x < m_stoneGridWidth && y >= 0 && y < m_stoneGridHeight) {
        int index = getStoneGridIndex(x, y);
        m_stoneGrid[index].occupied = false;
        m_stoneGrid[index].particleIndex = -1;
    }
}