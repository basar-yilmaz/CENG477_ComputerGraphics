#version 120

vec3 lightPos = vec3(-5, 5, -5);
uniform vec3 eyePos; // Make eyePos a uniform to be set from your application

vec3 I = vec3(0.8, 0.8, 0.1);

vec3 kd = vec3(0.1, 0.6, 0.3);
vec3 ka = vec3(0.1, 0.1, 0.1);
vec3 ks = vec3(0.1, 0.6, 0.3);

varying vec4 fragPos;
varying vec3 N;

void main(void)
{
    vec3 viewDir = normalize(eyePos - vec3(fragPos));
    vec3 lightDir = normalize(lightPos - vec3(fragPos));
    vec3 halfwayDir = normalize(lightDir + viewDir);

    float NdotL = dot(N, lightDir);
    float NdotH = dot(N, halfwayDir);

    vec3 diffuseColor = I * kd * max(0, NdotL);
    vec3 specularColor = I * ks * pow(max(0, NdotH), 20);

    gl_FragColor = vec4(diffuseColor + specularColor, 1);
}
