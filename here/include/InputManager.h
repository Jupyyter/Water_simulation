#pragma once
#include <SFML/Graphics.hpp>
#include <glm/glm.hpp>
#include "Config.h"
#include "FluidSystem.h"

class Camera;

class InputManager {
public:
    explicit InputManager(const Config& config);
    
    void handleEvent(const sf::Event& event, sf::RenderWindow& window, Camera* camera);
    void update(FluidSystem* fluidSystem, Camera* camera, const sf::RenderWindow& window);
    
private:
    void handleKeyPress(sf::Keyboard::Key key);
    
private:
    const Config& m_config;
    
    // Mouse interaction states
    bool m_isLeftMousePressed = false;
    bool m_isRightMousePressed = false; // For deletion
    bool m_isPanning = false;
    sf::Vector2i m_panStartPos{};
    glm::vec2 m_lastStonePos{-1.0f, -1.0f};
    
    // The current tool selected by the user
    MaterialType m_currentTool = MaterialType::WATER;
    
    // Radius for spawning and deleting
    float m_interactionRadius = 0.2f;
    int m_particlesPerClick = 5;
};