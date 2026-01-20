#version 330 core

in vec3 FragPos;
in vec3 Normal;

out vec4 FragColor;

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;
uniform vec3 objectColor;

void main()
{
    // Ambient
    float ambientStrength = 0.4;
    vec3 ambient = ambientStrength * lightColor;

    // Diffuse (smooth shading happens here via interpolated normals)
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);

    vec3 diffuse = 0.5*diff * lightColor;

    // Specular (Blinn-Phong)

    float specularStrength = 0.2;
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(norm, halfwayDir), 0.0), 32.0);
    vec3 specular = specularStrength * spec * lightColor;

    float fresnel = pow(1.0 - abs(dot(norm, viewDir)), 3.0);

    vec3 rimColor = vec3(0.5, 0.5, 0.5);



    vec3 result = (ambient + diffuse + specular) * objectColor + fresnel * rimColor;
//    vec3 result = (ambient + diffuse);
//    vec3 result = ambient;
//    vec3 result = vec3(1,1,1);
    FragColor = vec4(result, 1.0);
}