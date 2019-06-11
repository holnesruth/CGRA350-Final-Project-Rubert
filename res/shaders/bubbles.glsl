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
uniform sampler2D uNoise;
uniform samplerCube uCubeMap;

uniform vec2 uThickness;
uniform vec2 uLightEffects;

uniform float uTime;
uniform vec2 uResolution;

// viewspace data (this must match the output of the vertex shader)
in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} f_in;

// framebuffer output
out vec4 fb_color;

const vec3 DIRECTION = vec3(1, 1, 4);

/*============================ FLOW NOISE IMPLEMENTATION ============================*/
// With help from https://www.shadertoy.com/view/lslXRS

float time = uTime*0.1;

float hash21(vec2 n) { 
	return fract(sin(dot(n, vec2(12.9898, 4.1414)))*43758.5453); 
}

mat2 makem2(float theta) { 
	float c = cos(theta); 
	float s = sin(theta); 
	
	return mat2(c, -s, s, c); 
}

float noise(vec2 x) { 
	return texture(uNoise, x*0.1).x; 
}

vec2 gradient(vec2 p) {
	float ep = .09;
	float gradx = noise(vec2(p.x + ep, p.y)) - noise(vec2(p.x - ep, p.y));
	float grady = noise(vec2(p.x, p.y + ep)) - noise(vec2(p.x, p.y - ep));
	return vec2(gradx, grady);
}

float flow(vec2 p) {
	float z=2.;
	float rz = 0.;
	vec2 bp = p;
	for (float i = 1.0; i < 7.0; i++) {
		//primary flow speed
		p += time*0.6;
		
		//secondary flow speed (speed of the perceived flow)
		bp += time*1.9;
		
		//displacement field (try changing time multiplier)
		vec2 gr = gradient(i*p*0.34 + time*1.0);
		
		//rotation of the displacement field
		gr *= makem2(time*6.0 - (0.05*p.x + 0.03*p.y)*40.0);
		
		//displace the system
		p += gr*0.5;
		
		//add noise octave
		rz += (sin(noise(p)*7.0)*0.5 + 0.5)/z;
		
		//blend factor (blending displaced system with base system)
		//you could call this advection factor (.5 being low, .95 being high)
		p = mix(bp, p, 0.77);
		
		//intensity scaling
		z *= 1.4;
		//octave scaling
		p *= 2.0;
		bp *= 1.9;
	}
	return rz;	
}

/*============================ REFRACTION IMPLEMENTATION ============================*/

float snellsLaw(float theta_i, float n1, float n2) {
	return asin((n1/n2)*sin(theta_i));
}

float fresnel(float theta_i, float theta_t, float n1, float n2) {
	float num = n1*cos(theta_i) - n2*cos(theta_t);
	float den = n1*cos(theta_i) + n2*cos(theta_t);
	float Rs = pow(abs(num/den), 2);

	num = n1*cos(theta_t) - n2*cos(theta_i);
	den = n1*cos(theta_t) + n2*cos(theta_i);
	float Rp = pow(abs(num/den), 2);

	return 2.0*(Rs + Rp); //Scalar should be 0.5 but the colours aren't very vibrant
}

float calculateLightColor(float lambda, float thickness, float nAir, float nFilm, float theta_i, float intensity) {
	float theta_t = max(snellsLaw(theta_i, nAir, nFilm), 0.0);

	float d = 2 * M_PI * thickness * nFilm * cos(theta_t);
	float sind = sin(d/lambda);

	float fresnel = fresnel(theta_i, theta_t, nAir, nFilm);

	return 4.0 * intensity * fresnel * sind * sind;
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

	vec2 p = (f_in.position.xy / uResolution.xy) - 0.5;
	p.x *= uResolution.x/uResolution.y;
	p *= 3.0;
	float rz = flow(p);

	vec3 thickness = texture(uTexture, f_in.textureCoord).rgb;
	float t = (thickness.x + thickness.y + thickness.z)/3.0;
	float w = ((minThickness*(1.0 - t)) + (maxThickness*t)); 
	w = w/rz; //introduce noise

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