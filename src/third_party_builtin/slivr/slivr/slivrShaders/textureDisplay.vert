varying vec4 clip_Vertex;         

void main(){                      
	// get the texture coordinate  
	gl_TexCoord[0] = gl_TextureMatrix[0] * gl_MultiTexCoord0.xyzw; 

	// vertex as a vector - clip space coordinate                     
   clip_Vertex = vec4(gl_ModelViewProjectionMatrix * gl_Vertex);  

	// set the position of the vertex to unchanged!
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;//ftransform();
}