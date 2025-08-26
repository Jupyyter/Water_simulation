#include "InputManager.h"
#include "FluidSystem.h"
#include "Camera.h"
#include <iostream>
#include <algorithm>

InputManager::InputManager(const Config& config) 
    : m_config(config)
{
    std::cout << "InputManager initialized" << std::endl;
}

void InputManager::handleEvent(const sf::Event& event, sf::RenderWindow& window, Camera* camera) {
    if (const auto* mouseButtonPressed = event.getIf<sf::Event::MouseButtonPressed>()) {
        const sf::Vector2i mousePos = {mouseButtonPressed->position.x, mouseButtonPressed->position.y};

        switch (mouseButtonPressed->button) {
            case sf::Mouse::Button::Left:
                m_isLeftMousePressed = true;
                break;
            case sf::Mouse::Button::Right:
                m_isRightMousePressed = true;
                break;
            case sf::Mouse::Button::Middle:
                m_isPanning = true;
                m_panStartPos = mousePos;
                break;
            default:
                break;
        }
    }
    else if (const auto* mouseButtonReleased = event.getIf<sf::Event::MouseButtonReleased>()) {
        switch (mouseButtonReleased->button) {
            case sf::Mouse::Button::Left:
                m_isLeftMousePressed = false;
                break;
            case sf::Mouse::Button::Right:
                m_isRightMousePressed = false;
                break;
            case sf::Mouse::Button::Middle:
                m_isPanning = false;
                break;
            default:
                break;
        }
    }
    else if (const auto* mouseMoved = event.getIf<sf::Event::MouseMoved>()) {
        if (m_isPanning) {
            const sf::Vector2i newMousePos = {mouseMoved->position.x, mouseMoved->position.y};
            const sf::Vector2f lastWorldPos = camera->screenToWorld(m_panStartPos, window);
            const sf::Vector2f currentWorldPos = camera->screenToWorld(newMousePos, window);
            const sf::Vector2f delta = lastWorldPos - currentWorldPos;
            camera->pan(delta);
            m_panStartPos = newMousePos;
        }
    }
    else if (const auto* mouseWheelScrolled = event.getIf<sf::Event::MouseWheelScrolled>()) {
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::LControl) || sf::Keyboard::isKeyPressed(sf::Keyboard::Key::RControl)) {
            // Adjust interaction radius
            const float radiusFactor = (mouseWheelScrolled->delta > 0) ? 1.1f : 0.9f;
            m_interactionRadius = std::clamp(m_interactionRadius * radiusFactor, 0.05f, 1.0f);
            std::cout << "Interaction radius: " << m_interactionRadius << std::endl;
        } else {
            // Zoom camera
            const float zoomFactor = (mouseWheelScrolled->delta > 0) ? 1.f / 1.1f : 1.1f;
            camera->zoom(zoomFactor, mouseWheelScrolled->position, window);
        }
    }
    else if (const auto* keyPressed = event.getIf<sf::Event::KeyPressed>()) {
        handleKeyPress(keyPressed->code);
    }
}

void InputManager::handleKeyPress(sf::Keyboard::Key key) {
    switch (key) {
        case sf::Keyboard::Key::Num1: 
            m_currentTool = MaterialType::WATER; 
            std::cout << "Tool: Water\n"; 
            break;
        case sf::Keyboard::Key::Num2: 
            m_currentTool = MaterialType::SMOKE; 
            std::cout << "Tool: Smoke\n"; 
            break;
        case sf::Keyboard::Key::Num3: 
            m_currentTool = MaterialType::FIRE;  
            std::cout << "Tool: Fire\n"; 
            break;
        case sf::Keyboard::Key::Num4: 
            m_currentTool = MaterialType::SAND;  
            std::cout << "Tool: Sand\n"; 
            break;
        case sf::Keyboard::Key::Num5: 
            m_currentTool = MaterialType::OIL;   
            std::cout << "Tool: Oil\n"; 
            break;
        case sf::Keyboard::Key::Num6: 
            m_currentTool = MaterialType::STONE; 
            std::cout << "Tool: Stone\n"; 
            break;
        case sf::Keyboard::Key::Equal: 
            m_particlesPerClick = std::min(m_particlesPerClick + 10, 200);
            std::cout << "Particles per click: " << m_particlesPerClick << std::endl;
            break;
        case sf::Keyboard::Key::Hyphen:
            m_particlesPerClick = std::max(m_particlesPerClick - 10, 1);
            std::cout << "Particles per click: " << m_particlesPerClick << std::endl;
            break;
        default: 
            break;
    }
}

void InputManager::update(FluidSystem* fluidSystem, Camera* camera, const sf::RenderWindow& window) {
    if (!fluidSystem || !camera) return;

    const sf::Vector2f worldPosF = camera->screenToWorld(sf::Mouse::getPosition(window), window);
    const glm::vec2 worldPos = {worldPosF.x, worldPosF.y};

    // --- HANDLE LEFT CLICK (CREATION) ---
    if (m_isLeftMousePressed) {
        if (m_currentTool == MaterialType::STONE || m_currentTool == MaterialType::TERRAIN) {
            // For stone particles, calculate count based on area to prevent excessive creation
            const float area = m_interactionRadius * m_interactionRadius * 3.14159f;
            const float particleArea = m_config.particleRadius * m_config.particleRadius * 4.0f;
            int stoneCount = static_cast<int>(area / particleArea);
            stoneCount = std::clamp(stoneCount, 1, 50); // Reasonable limits
            
            fluidSystem->addParticles(worldPos, m_interactionRadius, stoneCount, MaterialType::STONE);
        } else {
            // Add particles for any other tool
            fluidSystem->addParticles(worldPos, m_interactionRadius, 
                                      m_particlesPerClick, m_currentTool);
        }
    }

    // --- HANDLE RIGHT CLICK (DELETION) ---
    if (m_isRightMousePressed) {
        // Delete particles for any tool (including stone)
        fluidSystem->removeParticlesInArea(worldPos, m_interactionRadius);
    }
}