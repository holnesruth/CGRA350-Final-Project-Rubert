#version 330 core

#define M_PI 3.1415926535897932384626433832795
#define ETA 1.4

#define RED 700
#define GREEN 560
#define BLUE 460

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uModelMatrix;

uniform sampler2D uTexture;
uniform samplerCube uCubeMap;

uniform vec2 uThickness;
uniform vec2 uLightEffects;

// viewspace data (this must match the output of the vertex shader)
in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} f_in;

// framebuffer output
out vec4 fb_color;

const vec3 DIRECTION = vec3(1, 1, 4);

float snellsLaw(float theta_i, float n1, float n2) {
	return asin((n1/n2) * sin(theta_i));
}

float fresnel(float theta_i, float theta_t, float n1, float n2) {
	float num = n1 * cos(theta_i) - n2 * cos(theta_t);
	float den = n1 * cos(theta_i) + n2 * cos(theta_t);
	float Rs = pow(abs(num / den), 2);

	num = n1 * cos(theta_t) - n2 * cos(theta_i);
	den = n1 * cos(theta_t) + n2 * cos(theta_i);
	float Rp = pow(abs(num / den), 2);

	return 5 * (Rs + Rp); //Scalar should be 0.5 but the colours aren't very vibrant
}

float calculateLightColor(float lambda, float thickness, float nAir, float nFilm, float theta_i, float intensity) {
	float theta_t = max(snellsLaw(theta_i, nAir, nFilm), 0.0);

	float d = 2 * M_PI * thickness * nFilm * cos(theta_t);
	float sind = sin(d / lambda);

	float fresnel = fresnel(theta_i, theta_t, nAir, nFilm);

	return 4.0 * intensity * fresnel * sind * sind;
}

float hash( float n ){
    return fract(sin(n)*43758.5453);
}

float noise( in vec3 x ) {
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f*f*(3.0-2.0*f);
    float n = p.x + p.y*57.0 + 113.0*p.z;
    return mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                   mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
               mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                   mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
}

vec3 noise3(vec3 x) {
	return vec3( noise(x+vec3(123.456,.567,.37)),
				 noise(x+vec3(.11,47.43,19.17)),
				 noise(x) );
}

void main() {
	float minThickness = uThickness.x;
	float maxThickness = uThickness.y;
	float intensity = uLightEffects.x;
	float transparency = uLightEffects.y;

	vec3 norm = f_in.normal;

	vec4 light = uProjectionMatrix * vec4(DIRECTION, 1.0);
	vec3 lightDir = normalize(vec3(light));

	vec3 viewDir = normalize(-f_in.position);
	vec3 reflectDir = reflect(-lightDir, norm);

	vec3 thickness = texture(uTexture, f_in.textureCoord).rgb + noise3(f_in.normal);
	float t = (thickness.x + thickness.y + thickness.z)/3.0;
	float w = ((minThickness * (1.0 - t)) + (maxThickness * t));

	float theta_i = max(dot(norm, lightDir), 0.0);

	float red = calculateLightColor(RED, w, 1, ETA, theta_i, intensity);
	float green = calculateLightColor(GREEN, w, 1, ETA, theta_i, intensity);
	float blue = calculateLightColor(BLUE, w, 1, ETA, theta_i, intensity);

	vec3 filmColour = vec3(red, green, blue);
	vec3 reflection = texture(uCubeMap, reflectDir).rgb;

	vec3 finalShading = ((1 - transparency) * filmColour) + (intensity * reflection);
	
	// output to the framebuffer
	fb_color = vec4(finalShading, transparency);
}