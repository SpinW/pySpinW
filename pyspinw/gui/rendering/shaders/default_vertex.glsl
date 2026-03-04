#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec3 aNormal;

uniform mat4 model;
uniform mat4 projectionView;


void main()
{
    gl_Position = projectionView * model * vec4(position, 1.0);
}