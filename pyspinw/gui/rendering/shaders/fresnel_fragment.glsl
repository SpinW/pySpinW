#version 410 core

in vec3 FragPos;
in vec3 Normal;

out vec4 FragColor;

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;
uniform vec3 objectColor;

void main()
{

    vec3 N = normalize(Normal);
    vec3 V = normalize(viewPos - FragPos);

    float fresnel = pow(1.0 - max(dot(N, V), 0.0), 4.0);

    vec3 color = vec3(0.2, 0.6, 1.0);
    float alpha = fresnel;

    FragColor = vec4(color * alpha, alpha);
}