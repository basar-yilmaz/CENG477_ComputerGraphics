#version 330 core

vec3 lightPos = vec3(-5, 5, -5);
vec3 eyePos = vec3(0, 0, 0);

uniform vec3 cubeColor;

// vec3 I = cubeColor;
// vec3 Iamb = vec3(0.8, 0.8, 0.8);

vec3 kd = vec3(0.8, 0.8, 0.3);
// vec3 ka = vec3(0.1, 0.1, 0.1);
vec3 ks = vec3(0.8, 0.8, 0.3);

in vec4 fragPos;
in vec3 N;

out vec4 outColor;

void main(void)
{
	vec3 L = normalize(lightPos - vec3(fragPos));
	vec3 V = normalize(eyePos - vec3(fragPos));
	vec3 H = normalize(L + V);
	float NdotL = dot(N, L);
	float NdotH = dot(N, H);

	vec3 diffuseColor = cubeColor * kd * max(0, NdotL);
	// vec3 ambientColor = Iamb * ka;
	vec3 specularColor = cubeColor * ks * pow(max(0, NdotH), 20);

    outColor = vec4(diffuseColor + specularColor, 1);
}
