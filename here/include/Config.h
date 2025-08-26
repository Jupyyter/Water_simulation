#pragma once
#include <box2d/box2d.h>

struct Config {
    // Box2D settings
    float box2dGravity = -9.81f;
    int velocityIterations = 8;
    int positionIterations = 3;
    float box2dTimeStep = 1.0f / 60.0f;
    
    // Particle collision settings
    bool useBox2DForParticles = true;
    float particleRestitution = 0.1f;
    float particleFriction = 0.3f;
    float particleDensityB2D = 1.0f;
    
    // Existing settings
    bool separateParticles = true;
    bool compensateDrift = true;
    int windowWidth = 1200;
    int windowHeight = 800;
    const char* windowTitle = "2D Physics Simulator";
    int targetFPS = 60;
    
    float physicsTimeStep = 1.0f / 60.0f;
    float gravity = -9.81f;
    int maxParticles = 10000;
    
    float simWidth = 12.0f;
    float simHeight = 8.0f;
    
    float particleRadius = 0.03f;
    float fluidDensity = 1000.0f;
    float flipRatio = 0.8f; // Reduced from 0.9f for better stability
    int pressureIterations = 60; // Reduced from 80 for better performance
    int particleSeparationIterations = 3; // Reduced from 4 for better performance
    
    bool showDebugUI = true;
    bool showParticles = true;
    bool showGrid = false;
    bool enableVSync = true;
    
    // Cached scale value for performance
    mutable float m_cachedScale = -1.0f;
    
    float getScale() const {
        if (m_cachedScale < 0.0f) {
            m_cachedScale = windowHeight / simHeight;
        }
        return m_cachedScale;
    }
};