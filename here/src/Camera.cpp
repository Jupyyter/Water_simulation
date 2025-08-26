#include "Camera.h"
#include "FluidSystem.h"
#include <algorithm>
#include <iostream>
#include <cmath>

Camera::Camera(const Config& config) 
    : m_config(config)
{
    m_view.setSize({config.simWidth, -config.simHeight});
    m_view.setCenter({config.simWidth / 2.0f, config.simHeight / 2.0f});

    float visualRadius = config.particleRadius * 0.8f;
    m_particleShape.setRadius(visualRadius);
    m_particleShape.setOrigin({visualRadius, visualRadius});
    
    float cellSize = (config.simHeight / 100.0f);
    m_cellShape.setSize({cellSize, cellSize});
    m_cellShape.setOrigin({cellSize * 0.5f, cellSize * 0.5f});
    
    m_obstacleShape.setFillColor(sf::Color::Red);
    m_terrainShape.setFillColor(sf::Color(128, 128, 128));

    if (m_font.openFromFile("assets/fonts/arial.ttf")) {
        m_fontLoaded = true;
        m_text.setFont(m_font);
        m_text.setCharacterSize(12);
        m_text.setFillColor(sf::Color::White);
    }
    
    std::cout << "Camera initialized" << std::endl;
}

void Camera::applyTo(sf::RenderWindow& window) const {
    window.setView(m_view);
}

sf::Vector2f Camera::screenToWorld(const sf::Vector2i& screenPos, const sf::RenderWindow& window) const {
    return window.mapPixelToCoords(screenPos, m_view);
}

void Camera::zoom(float factor, const sf::Vector2i& mousePixelPos, sf::RenderWindow& window) {
    const sf::Vector2f beforeCoord = window.mapPixelToCoords(mousePixelPos, m_view);
    m_view.zoom(factor);
    const sf::Vector2f afterCoord = window.mapPixelToCoords(mousePixelPos, m_view);
    const sf::Vector2f offset = beforeCoord - afterCoord;
    m_view.move(offset);
}

void Camera::pan(const sf::Vector2f& delta) {
    m_view.move(delta);
}

void Camera::render(sf::RenderWindow& window, const FluidSystem& fluidSystem, const Config& config) {
    
    if (config.showGrid) {
        renderGrid(window, fluidSystem);
    }
    
    if (config.showParticles) {
        renderParticles(window, fluidSystem);
    }
}

void Camera::renderParticles(sf::RenderWindow& window, const FluidSystem& fluidSystem) {
    const auto& particles = fluidSystem.getParticles();

    for (const auto& particle : particles) {
        if (!particle.active) continue;
        
        m_particleShape.setPosition({particle.position.x, particle.position.y});
        
        sf::Color color = vec3ToColor(particle.color);
        
        if (particle.material == MaterialType::SMOKE || particle.material == MaterialType::FIRE) {
            color.a = static_cast<std::uint8_t>(255 * particle.life);
        }
        
        m_particleShape.setFillColor(color);
        window.draw(m_particleShape);
    }
}

void Camera::renderGrid(sf::RenderWindow& window, const FluidSystem& fluidSystem) {
    const auto& cellTypes = fluidSystem.getCellTypes();
    const auto& gridData = fluidSystem.getGrid();
    
    int gridWidth = fluidSystem.getGridWidth();
    int gridHeight = fluidSystem.getGridHeight();
    float cellSize = fluidSystem.getCellSize();
    
    m_cellShape.setSize({cellSize, cellSize});
    m_cellShape.setOrigin({cellSize * 0.5f, cellSize * 0.5f});
    
    for (int i = 0; i < gridWidth; ++i) {
        for (int j = 0; j < gridHeight; ++j) {
            int index = i * gridHeight + j;
            CellType cellType = cellTypes[index];
            
            if (cellType == CellType::AIR) continue;
            
            glm::vec2 worldPos = glm::vec2((i + 0.5f) * cellSize, (j + 0.5f) * cellSize);
            m_cellShape.setPosition({worldPos.x, worldPos.y});
            
            sf::Color cellColor;
            switch (cellType) {
                case CellType::FLUID:
                    cellColor = densityToColor(gridData[index]);
                    cellColor.a = 128;
                    break;
                case CellType::SOLID:
                    cellColor = sf::Color(80, 80, 80, 200);
                    break;
                case CellType::AIR:
                default:
                    cellColor = sf::Color::Transparent;
                    break;
            }
            
            m_cellShape.setFillColor(cellColor);
            m_cellShape.setOutlineThickness(0.01f);
            m_cellShape.setOutlineColor(sf::Color(50, 50, 50, 100));
            
            window.draw(m_cellShape);
        }
    }
}


sf::Color Camera::vec3ToColor(const glm::vec3& color, float alpha) const {
    return sf::Color(
        static_cast<std::uint8_t>(glm::clamp(color.r * 255.0f, 0.0f, 255.0f)),
        static_cast<std::uint8_t>(glm::clamp(color.g * 255.0f, 0.0f, 255.0f)),
        static_cast<std::uint8_t>(glm::clamp(color.b * 255.0f, 0.0f, 255.0f)),
        static_cast<std::uint8_t>(glm::clamp(alpha * 255.0f, 0.0f, 255.0f))
    );
}

sf::Color Camera::densityToColor(float density, float maxDensity) const {
    float normalizedDensity = glm::clamp(density / maxDensity, 0.0f, 1.0f);
    
    std::uint8_t red = static_cast<std::uint8_t>(normalizedDensity * 255);
    std::uint8_t blue = static_cast<std::uint8_t>((1.0f - normalizedDensity) * 255);
    std::uint8_t green = static_cast<std::uint8_t>(std::abs(0.5f - normalizedDensity) * 2.0f * 255);
    
    return sf::Color(red, green, blue, 128);
}