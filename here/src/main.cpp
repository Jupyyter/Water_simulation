#include <SFML/Graphics.hpp>
#include <iostream>
#include <memory>

#include "PhysicsSimulator.h"
#include "Config.h"

int main() {
    try {
        // Initialize configuration
        Config config;
        config.windowWidth = 1200;
        config.windowHeight = 800;
        config.windowTitle = "2D Physics Simulator";
        config.targetFPS = 60;
        config.physicsTimeStep = 1.0f / 120.0f; // Physics runs at 120Hz
        
        // Create and run the simulator
        PhysicsSimulator simulator(config);
        simulator.run();
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }
    
    return 0;
}