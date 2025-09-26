#pragma once

#include <vector>
#include <memory>
#include <glm/glm.hpp>
#include <box2d/box2d.h>
#include "Config.h"

enum class MaterialType {
    WATER = 0,
    SMOKE,
    FIRE,
    SAND,
    OIL,
    STONE,
    TERRAIN // Keep for backwards compatibility, but will be treated as STONE
};

struct Particle {
    glm::vec2 position{0.0f};
    glm::vec2 velocity{0.0f};
    glm::vec3 color{0.0f, 0.4f, 1.0f};
    
    // Physical properties
    float mass = 1.0f;
    float density = 1.0f;
    float pressure = 0.0f;
    float temperature = 20.0f;

    // Simulation state
    MaterialType material = MaterialType::WATER;
    float life = 1.0f;
    bool active = true;
    bool isStatic = false;
};

enum class CellType { FLUID = 0, AIR, SOLID };

class FluidSystem {
public:
    explicit FluidSystem(const Config& config);
    ~FluidSystem();
    void enforceGridBoundary();
    void update(float deltaTime);
    void reset();
    
    void addParticles(const glm::vec2& position, float radius, int count, MaterialType material = MaterialType::WATER);
    void removeParticlesInArea(const glm::vec2& position, float radius);
    
    void addStaticTerrain(const glm::vec2& position, float radius);
    void removeTerrainInArea(const glm::vec2& position, float radius);
    
    // Getters
    const std::vector<Particle>& getParticles() const { return m_particles; }
    const std::vector<float>& getGrid() const { return m_grid; }
    const std::vector<CellType>& getCellTypes() const { return m_cellTypes; }
    int getParticleCount() const;
    int getGridWidth() const { return m_gridWidth; }
    int getGridHeight() const { return m_gridHeight; }
    float getCellSize() const { return m_cellSize; }

private:
    void applyGravity(float dt);
    void handleBoundaryCollisions();
    void separateParticles(int iterations);
    void transferVelocitiesToGrid();
    void solveIncompressibility(int iterations, float dt);
    void transferVelocitiesFromGrid();
    void updateParticleColors();
    void buildSpatialHash();
    void buildStoneHash(); // Separate hash just for stone particles
    void handleStoneCollisions();
    void updateSolidCellsFromStones();

    inline int getGridIndex(int x, int y) const { return x * m_gridHeight + y; }
    glm::vec2 getGridPosition(int x, int y) const;
    glm::ivec2 worldToGrid(const glm::vec2& worldPos) const;
    glm::vec2 gridToWorld(const glm::ivec2& gridPos) const;
    
    // Stone grid functions
    glm::ivec2 worldToStoneGrid(const glm::vec2& worldPos) const;
    glm::vec2 stoneGridToWorld(const glm::ivec2& gridPos) const;
    inline int getStoneGridIndex(int x, int y) const { return x * m_stoneGridHeight + y; }
    inline bool isStoneGridCellOccupied(int x, int y) const;
    inline void addStoneToGrid(int x, int y, int particleIndex);
    inline void removeStoneFromGrid(int x, int y);
    
    float interpolateVelocityComponent(float x, float y, const std::vector<float>& field, float offsetX, float offsetY) const;
    void transferVelocityComponent(float x, float y, float vel, float mass, std::vector<float>& field, std::vector<float>& weights, float offsetX, float offsetY);
    
    void updateSmokeParticles(float dt);
    void updateFireParticles(float dt);
    void cleanupInactiveParticles();
    
private:
    const Config& m_config;
    
    std::vector<Particle> m_particles;
    std::vector<int> m_activeParticles;
    std::vector<int> m_stoneParticles;
    
    // Stone grid system for non-overlapping placement
    struct StoneGridCell {
        bool occupied = false;
        int particleIndex = -1; // Index of the stone particle in this cell
    };
    std::vector<StoneGridCell> m_stoneGrid;
    int m_stoneGridWidth, m_stoneGridHeight;
    float m_stoneGridCellSize; // Size of each stone grid cell (should be ~particle diameter)
    
    int m_gridWidth, m_gridHeight;
    float m_cellSize;
    float m_invCellSize;
    
    std::vector<float> m_velocityU, m_velocityV;
    std::vector<float> m_prevVelocityU, m_prevVelocityV;
    std::vector<float> m_pressure;
    std::vector<float> m_divergence;
    std::vector<float> m_particleDensity;
    std::vector<CellType> m_cellTypes;
    std::vector<float> m_grid;
    
    // Grids for mass-weighted velocity transfer
    std::vector<float> m_massU, m_massV;

    // Particle spatial hash (only for fluid particles)
    struct SpatialCell { std::vector<int> particles; };
    std::vector<SpatialCell> m_spatialHash;
    
    // Separate stone spatial hash (updated less frequently)
    std::vector<SpatialCell> m_stoneHash;
    bool m_stoneHashNeedsUpdate = true; // Flag to update stone hash only when needed
    
    int m_hashGridWidth, m_hashGridHeight;
    float m_hashCellSize;
    float m_invHashCellSize; // Precomputed inverse for optimization

    float m_particleRestDensity = 0.0f;
    float m_overRelaxation = 1.9f;
    
    struct MaterialProperties {
        float mass;
        glm::vec3 color;
    };
    std::vector<MaterialProperties> m_materialProperties;
    
    // Cached values for optimization
    float m_minDist;
    float m_minDist2;
    float m_extraSeparation;
    
    void initializeMaterialProperties();
    void setupInitialFluid();
};