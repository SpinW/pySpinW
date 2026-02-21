#version 330 core

/* Render the ID an int buffer */

uniform uint objectID;

layout (location = 0) out uint FragID;

void main()
{
    FragID = objectID;
}