#pragma once

#include <SFML/Graphics.hpp>
#include <memory>

#include "Config.h"
#include "FluidSystem.h"
#include "Camera.h"
#include "InputManager.h"

class PhysicsSimulator {
public:
    explicit PhysicsSimulator(const Config& config);
    ~PhysicsSimulator() = default;
    
    void run();
    
private:
    void handleEvents();
    void update(float deltaTime);
    void render();
    void updateDebugUI();
    
    // Fixed timestep physics update
    void physicsUpdate(float fixedDeltaTime);
    
private:
    Config m_config;
    sf::RenderWindow m_window;
    sf::Clock m_deltaClock;
    
    // Core systems
    std::unique_ptr<FluidSystem> m_fluidSystem;
    std::unique_ptr<Camera> m_camera;
    std::unique_ptr<InputManager> m_inputManager;
    
    // Timing
    float m_physicsAccumulator = 0.0f;
    bool m_isPaused = false;
    
    // Debug
    sf::Font m_debugFont;
    sf::Text m_debugText{m_debugFont};
    int m_frameCount = 0;
    float m_fpsTimer = 0.0f;
    float m_currentFPS = 0.0f;
};