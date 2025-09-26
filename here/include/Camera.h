#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include <glm/glm.hpp>
#include "Config.h"

class FluidSystem;
struct Particle;

class Camera {
public:
    explicit Camera(const Config& config);
    ~Camera() = default;
    
    void render(sf::RenderWindow& window, const FluidSystem& fluidSystem, const Config& config);

    // Camera manipulation
    void zoom(float factor, const sf::Vector2i& mousePixelPos, sf::RenderWindow& window);
    void pan(const sf::Vector2f& delta);
    void applyTo(sf::RenderWindow& window) const;
    sf::Vector2f screenToWorld(const sf::Vector2i& screenPos, const sf::RenderWindow& window) const;
    const sf::View& getView() const { return m_view; }

private:
    void renderParticles(sf::RenderWindow& window, const FluidSystem& fluidSystem);
    void renderGrid(sf::RenderWindow& window, const FluidSystem& fluidSystem);

    sf::Color vec3ToColor(const glm::vec3& color, float alpha = 1.0f) const;
    sf::Color densityToColor(float density, float maxDensity = 2.0f) const;

private:
    const Config& m_config;
    sf::View m_view;
    
    sf::CircleShape m_particleShape;
    
    std::vector<sf::Vertex> m_gridVertices;
    sf::RectangleShape m_cellShape;
    
    sf::CircleShape m_obstacleShape;

    sf::CircleShape m_terrainShape;
    
    // Debug rendering
    sf::Font m_font;
    sf::Text m_text{m_font};
    bool m_fontLoaded = false;
};