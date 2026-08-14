#version 410 core

in vec3 FragPos;
in vec3 Normal;

out vec4 FragColor;

uniform vec3 viewPos;
uniform vec3 objectColor;

void main()
{

    vec3 N = normalize(Normal);
    vec3 V = normalize(viewPos - FragPos);

    float alpha = pow(1.0 - max(dot(N, V), 0.0), 2.0);

    FragColor = vec4(objectColor, alpha);
}