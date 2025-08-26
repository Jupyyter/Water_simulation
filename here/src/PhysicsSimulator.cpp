#include "PhysicsSimulator.h"
#include <iostream>
#include <algorithm>

PhysicsSimulator::PhysicsSimulator(const Config& config) 
    : m_config(config)
    , m_window(sf::VideoMode({static_cast<unsigned int>(config.windowWidth), static_cast<unsigned int>(config.windowHeight)}), config.windowTitle)
{
    if (m_config.enableVSync) {
        m_window.setVerticalSyncEnabled(true);
    } else {
        m_window.setFramerateLimit(m_config.targetFPS);
    }
    
    // Initialize systems
    m_fluidSystem = std::make_unique<FluidSystem>(m_config);
    m_camera = std::make_unique<Camera>(m_config);
    m_inputManager = std::make_unique<InputManager>(m_config);
    
    // Setup debug text
    if (!m_debugFont.openFromFile("assets/fonts/ARIAL.TTF")) {
        std::cout << "Warning: Could not load debug font, using default font\n";
    }
    m_debugText.setFont(m_debugFont);
    m_debugText.setCharacterSize(16);
    m_debugText.setFillColor(sf::Color::White);
    m_debugText.setPosition({10.f, 10.f});
    
    std::cout << "Physics Simulator initialized\n";
}

void PhysicsSimulator::run() {
    std::cout << "Starting simulation...\n";
    
    while (m_window.isOpen()) {
        float deltaTime = m_deltaClock.restart().asSeconds();
        
        // Cap delta time to prevent spiral of death
        deltaTime = std::min(deltaTime, 0.05f);
        
        handleEvents();
        update(deltaTime);
        render();
    }
}

void PhysicsSimulator::handleEvents() {
    while (auto event = m_window.pollEvent()) {
        if (event->is<sf::Event::Closed>()) {
            m_window.close();
        }
        
        // Handle input, now passing camera for view controls
        m_inputManager->handleEvent(*event, m_window, m_camera.get());
        
        // Global controls
        if (const auto* keyPressed = event->getIf<sf::Event::KeyPressed>()) {
            switch (keyPressed->code) {
                case sf::Keyboard::Key::Space:
                    m_isPaused = !m_isPaused;
                    std::cout << (m_isPaused ? "Paused" : "Resumed") << "\n";
                    break;
                case sf::Keyboard::Key::R:
                    m_fluidSystem->reset();
                    std::cout << "Simulation reset\n";
                    break;
                case sf::Keyboard::Key::G:
                    m_config.showGrid = !m_config.showGrid;
                    break;
                case sf::Keyboard::Key::P:
                    m_config.showParticles = !m_config.showParticles;
                    break;
                case sf::Keyboard::Key::Escape:
                    m_window.close();
                    break;
            }
        }
    }
    
    // Handle continuous input, now passing camera and window
    m_inputManager->update(m_fluidSystem.get(), m_camera.get(), m_window);
}

void PhysicsSimulator::update(float deltaTime) {
    if (m_isPaused) return;
    
    // Fixed timestep physics
    m_physicsAccumulator += deltaTime;
    
    while (m_physicsAccumulator >= m_config.physicsTimeStep) {
        physicsUpdate(m_config.physicsTimeStep);
        m_physicsAccumulator -= m_config.physicsTimeStep;
    }
    
    // Update FPS counter
    m_frameCount++;
    m_fpsTimer += deltaTime;
    if (m_fpsTimer >= 1.0f) {
        m_currentFPS = m_frameCount / m_fpsTimer;
        m_frameCount = 0;
        m_fpsTimer = 0.0f;
    }
}

void PhysicsSimulator::physicsUpdate(float fixedDeltaTime) {
    m_fluidSystem->update(fixedDeltaTime);
}

void PhysicsSimulator::render() {
    m_window.clear(sf::Color::Black);
    
    // Apply the camera's view for rendering the simulation world
    m_camera->applyTo(m_window);

    // Render systems using the camera
    m_camera->render(m_window, *m_fluidSystem, m_config);
    
    // Switch back to the default view to render UI elements
    m_window.setView(m_window.getDefaultView());
    
    // Render debug info
    if (m_config.showDebugUI) {
        updateDebugUI();
        m_window.draw(m_debugText);
    }
    
    m_window.display();
}

void PhysicsSimulator::updateDebugUI() {
    std::string debugInfo = 
        "FPS: " + std::to_string(static_cast<int>(m_currentFPS)) + "\n" +
        "Particles: " + std::to_string(m_fluidSystem->getParticleCount()) + "\n" +
        "Status: " + (m_isPaused ? "PAUSED" : "RUNNING") + "\n\n" +
        "Controls:\n" +
        "LMB - Add Material/Terrain\n" +
        "RMB - Delete Material/Terrain\n" +
        "MMB Drag - Pan camera\n" +
        "Mouse Wheel - Zoom camera\n" +
        "Ctrl + Wheel - Adjust spawn/delete radius\n" +
        "1-5 - Select Material | 6 - Terrain Tool\n" +
        "SPACE - Pause | R - Reset | ESC - Exit";
    
    m_debugText.setString(debugInfo);
}