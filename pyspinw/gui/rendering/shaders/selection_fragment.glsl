#version 330 core
out vec4 FragColor;

uniform vec3 selectionColor;

void main()
{
    FragColor = vec4(selectionColor, 1.0);
}