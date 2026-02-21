#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

uniform mat4 model;
uniform mat4 projectionView;

out vec3 FragPos;
out vec3 Normal;

void main()
{
    FragPos = vec3(model * vec4(aPos, 1.0)); // Position in world
    Normal = mat3(transpose(inverse(model))) * aNormal; // Normal

    gl_Position = projectionView * vec4(FragPos, 1.0);
}