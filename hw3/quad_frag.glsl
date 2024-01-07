#version 330 core

in vec4 fragPos; // The position in world space

uniform float scale; // Controls the size of the checkerboard squares
uniform vec3 color1; // First color (e.g., black)
uniform vec3 color2; // Second color (e.g., white)

out vec4 outColor;

void main() {
    // Apply the checkerboard pattern
    vec2 checker = floor(fragPos.xz * scale); // Modify the scale to control the size of the squares
    float checkerPattern = mod(checker.x + checker.y, 2.0); // Checkerboard pattern

    // Choose the color based on the checker pattern
    outColor = mix(vec4(color1, 1.0), vec4(color2, 1.0), step(0.5, checkerPattern));
}