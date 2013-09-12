varying vec4 clip_Vertex;      
varying vec2 transformed_c;    
 
void main(){                   
	// get the texture coordinate                     
	gl_TexCoord[0] = gl_MultiTexCoord0.xyzw;        
 
 	//Equation 12 in occlusion shading paper          
	vec4 sample_position = gl_ModelViewProjectionMatrix*vec4(1.0, 1.0, gl_Vertex.z, 1.0);  
	transformed_c = sample_position.xy/sample_position.w;		    
 
 	// vertex as a vector - clip space coordinate     
    clip_Vertex = vec4(gl_ModelViewProjectionMatrix * gl_Vertex);  
 
 	// set the position of the vertex to unchanged!   
	gl_Position = ftransform();                     
}