#version 330 core
layout (location = 0) in vec3 position;

uniform mat4 projectionView;

void main()
{
    gl_Position = projectionView * vec4(position, 1.0);
}