//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2008 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : ShaderProgramARB.cc
//    Author : Milan Ikits
//    Date   : Wed Jul  7 23:21:33 2004

#include <slivr/gldefs.h>
#include <slivr/ShaderProgramARB.h>
#include <slivr/Utils.h>
#include <slivr/slivr_debug.h>
#include <cstring>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>

using std::cerr;
using std::endl;
using std::string;

namespace SLIVR {

const char *TEX_DISP_VERT =
    "void main(){                       \n " \
    " // get the texture coordinate   \n " \
    " gl_TexCoord[0] = gl_TextureMatrix[0] * gl_MultiTexCoord0.xyzw;  \n " \
    " \n " \
    " gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;//ftransform(); \n " \
    "} \n ";

bool ShaderProgramARB::init_ = false;
bool ShaderProgramARB::supported_ = false;
bool ShaderProgramARB::non_2_textures_ = false;
int ShaderProgramARB::max_texture_size_1_ = 64;
int ShaderProgramARB::max_texture_size_4_ = 64;

const char *tex_strings[] = {"tex0", "tex1", "tex2", "tex3",
                                 "tex4", "tex5", "tex6", "tex7",
                                 "tex8", "tex9", "tex10", "tex11",
                                 "tex12", "tex13", "tex14", "tex15"};

const char *loc_strings[] = {"loc0", "loc1", "loc2", "loc3",
                                 "loc4", "loc5", "loc6", "loc7",
                                 "loc8", "loc9", "loc10", "loc11",
                                 "loc12", "loc13", "loc14", "loc15"};

//static Mutex ShaderProgramARB_init_Mutex("ShaderProgramARB Init Lock");  

ShaderProgramARB::ShaderProgramARB(const string& program) :
  type_(0), id_(0), program_(program)
{
}

ShaderProgramARB::ShaderProgramARB(const string& programFrag, const string& programVert) :
  type_(0), id_(0), programFrag_(programFrag), programVert_(programVert)
{
}

ShaderProgramARB::~ShaderProgramARB ()
{}

bool
ShaderProgramARB::valid()
{
  return shaders_supported() ? glIsProgram(id_) : false;
}

bool
ShaderProgramARB::init_shaders_supported(std::string& error)
{
  if (!init_)
  {
    const GLenum err = glewInit();
    if(err != GLEW_OK) {
        fprintf(stderr, "GLEW initialization failed: %s\n",
                glewGetErrorString(err));
        return false;
    }
    glewExperimental = GL_TRUE;

    CHECK_OPENGL_ERROR();

    if( glTexImage3D == NULL ) 
    {
      // Removed exit() statement, a library should never call exit()
      // Report an error and let the user of the slivr library decide what
      // to do. 
      error = "The OpenGL library's version is lower than 2.0.\n Please update your OpenGL graphic card drivers.\n";
#if DEBUG
      std::cerr << "(ShaderProgramARB::init_shaders_supported: first glTexImage3D check fails) " << error;
#endif
      return (false);
    }

    // TODO: Frame buffer objects are technically only required for 2D
    // transfer function generation and not 1D colormaps.  However the code
    // just blows up if you try to do a 2D colormap and they aren't supported.
    supported_ = GLEW_ARB_shading_language_100 && GLEW_EXT_framebuffer_object;

    // Blacklist Intel chipsets everywhere.  The drivers claim that they
    // support shaders but crash when you try to use them.  This
    // covers the Intel integrated chipsets in most laptops.
    const GLubyte* glRendererString = glGetString(GL_RENDERER);
    if (strncmp((const char *)glRendererString, "Intel", 5) == 0)
    {
      supported_ = false;
    }

#ifndef __sgi
    // Supported used to mean just shader programs.  However the
    // non-shader render path has been removed and so there is
    // currently no reason to check for texture sizes if they are not
    // supported.
    if (supported_)
    {
      int i;

      // Clear the OpenGL errors before checking for proxy textures.
      CHECK_OPENGL_ERROR();
      for (i = 64; i <= GL_MAX_TEXTURE_SIZE; i*=2)
      {
        GLint width = 0;
        if( glTexImage3D == NULL ) {
          const GLubyte* glVersionString;
          glVersionString = glGetString(GL_VERSION);

          // Removed exit() statement, a library should never call exit()
          // Report an error and let the user of the slivr library decide what
          // to do. 
          error = "The OpenGL library's version is lower than 2.0.\n Please update your OpenGL graphic card drivers.\n";
#if DEBUG
      std::cerr << "(ShaderProgramARB::init_shaders_supported: second glTexImage3D check fails) " << error;
#endif
          return (false);
        }
        glTexImage3D(GL_PROXY_TEXTURE_3D, 0, GL_LUMINANCE, i, i, i, 0,
                     GL_LUMINANCE, GL_UNSIGNED_BYTE, NULL);

        if (glGetError() == GL_NO_ERROR)
        {
          glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0,
                                   GL_TEXTURE_WIDTH, &width);
        }

        if (width == 0)
        {
          i /= 2;
          break;
        }
      }
      max_texture_size_1_ = Clamp(i, 64, i/2);


      // Clear the OpenGL errors before checking for proxy textures.
      CHECK_OPENGL_ERROR();
      for (i = 64; i <= GL_MAX_TEXTURE_SIZE; i*=2)
      {
        GLint width = 0;
        glTexImage3D(GL_PROXY_TEXTURE_3D, 0, GL_RGBA, i, i, i, 0,
                     GL_RGBA, GL_UNSIGNED_BYTE, NULL);
        if (glGetError() == GL_NO_ERROR)
        {
          glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0,
                                   GL_TEXTURE_WIDTH, &width);
        }

        if (width == 0)
        {
          i /= 2;
          break;
        }
      }
      max_texture_size_4_ = Clamp(i, 64, i/2);
    }
#endif // !sgi

    non_2_textures_ = GLEW_ARB_texture_non_power_of_two;
    init_ = true;
  }
  return (true);
}
  

bool
ShaderProgramARB::shaders_supported()
{
  return supported_;
}

bool
ShaderProgramARB::initialized()
{
  return init_;
}


int
ShaderProgramARB::max_texture_size_1()
{
  return max_texture_size_1_;
}

int
ShaderProgramARB::max_texture_size_4()
{
  return max_texture_size_4_;
}

bool
ShaderProgramARB::texture_non_power_of_two()
{
  return non_2_textures_;
}

bool
ShaderProgramARB::create(std::string& error)
{
  CHECK_OPENGL_ERROR();

  if (shaders_supported())
  {
    GLuint shader;
    switch(type_) 
    {
      case GL_VERTEX_PROGRAM_ARB:
        shader = glCreateShader(GL_VERTEX_SHADER);
        break;
      case GL_FRAGMENT_PROGRAM_ARB:
        shader = glCreateShader(GL_FRAGMENT_SHADER);
        break;
      default:
        error = "Error assigning shader type.";
        return (false);
    }
    
    if (shader == 0) 
    {
      error = "Error creating shader handle.";
      return (false);
    }

    //    cerr << program_ << endl;

    // set the source code and compile the shader
    const char *source[1];
    source[0] = program_.c_str();
    glShaderSource(shader, 1, source, NULL);
    glCompileShader(shader);
    
    // check the compilation of the shader
    GLint shader_status[1];
    glGetShaderiv(shader, GL_COMPILE_STATUS, shader_status);
    if (shader_status[0] == GL_FALSE) 
    {
      std::ostringstream oss;
      oss << "Error compiling shader.\n";
      
      char shader_log[1000];
      GLint shader_length[1];
      glGetShaderInfoLog(shader, 1000, shader_length, shader_log);
      
      oss << shader_log << endl;

      int line = 1;
      oss << std::setw(3) << line++ << ":";
      for (unsigned int i = 0; i < program_.length(); i++) 
      {
        if (program_[i] == '\n') 
        {
          oss << std::setw(3) << endl << line++ << ":";
        }
        else 
        {
          oss << program_[i];
        }
      }
      
      // A library should never call exit() 
      error = oss.str();
      return (false);
    }
    
    // create the GLSL program and attach the shader
    id_ = glCreateProgram();
    if (id_ == 0) 
    {
      error = "Error creating GLSL program";
      return (false);
    }
    glAttachShader(id_, shader);
    return (true);
  }
  
  CHECK_OPENGL_ERROR();
  return (false);
}

bool
ShaderProgramARB::createBoth(std::string& error)
{
  CHECK_OPENGL_ERROR();
  if (shaders_supported())
  {
    GLuint shaderVert, shaderFrag;
    shaderVert = glCreateShader(GL_VERTEX_SHADER);
    shaderFrag = glCreateShader(GL_FRAGMENT_SHADER);
      
    if (shaderVert == 0 || shaderFrag == 0) 
    {
      error = "Error creating shader handle.";
      return (false);
    }

    // create the GLSL program and attach the shader
    // set the source code and compile the shader
    const char *sourceVert[1], *sourceFrag[1];
    sourceVert[0] = programVert_.c_str();
    glShaderSource(shaderVert, 1, sourceVert, NULL);
    glCompileShader(shaderVert);
    
    // check the compilation of the shader
    GLint shader_status[1];
    glGetShaderiv(shaderVert, GL_COMPILE_STATUS, shader_status);
    if (shader_status[0] == GL_FALSE) 
    {
      std::ostringstream oss;
      oss << "Error compiling vertex shader.\n";
      
      char shader_log[1000];
      GLint shader_length[1];
      glGetShaderInfoLog(shaderVert, 1000, shader_length, shader_log);
      
      oss << shader_log << endl;

      int line = 1;
      oss << std::setw(3) << line++ << ":";
      for (unsigned int i = 0; i < program_.length(); i++) 
      {
        if (program_[i] == '\n') 
        {
          oss << std::setw(3) << endl << line++ << ":";
        }
        else 
        {
          oss << program_[i];
        }
      }
      
      // A library should never call exit() 
      error = oss.str();
      return (false);
    }

    // set the source code and compile the shader
    sourceFrag[0] = programFrag_.c_str();
    glShaderSource(shaderFrag, 1, sourceFrag, NULL);
    glCompileShader(shaderFrag);
    
    // check the compilation of the shader
    glGetShaderiv(shaderFrag, GL_COMPILE_STATUS, shader_status);
    if (shader_status[0] == GL_FALSE) 
    {
      std::ostringstream oss;
      oss << "Error compiling fragment shader.\n";
      
      char shader_log[1000];
      GLint shader_length[1];
      glGetShaderInfoLog(shaderFrag, 1000, shader_length, shader_log);
      
      oss << shader_log << endl;

      int line = 1;
      oss << std::setw(3) << line++ << ":";
      for (unsigned int i = 0; i < program_.length(); i++) 
      {
        if (program_[i] == '\n') 
        {
          oss << std::setw(3) << endl << line++ << ":";
        }
        else 
        {
          oss << program_[i];
        }
      }
      
      // A library should never call exit() 
      error = oss.str();
      return (false);
    }

    id_ = glCreateProgram();
    if (id_ == 0) 
    {
      error = "Error creating GLSL program";
      std::cout << " Error creating GLSL program " << std::endl;
      return (false);
    }

    glAttachShader(id_, shaderVert);
    glAttachShader(id_, shaderFrag);

    return (true);
  }

  CHECK_OPENGL_ERROR();
  return (false);
}


void
ShaderProgramARB::destroy ()
{
  if (shaders_supported())
  {
    glDeleteProgram(id_);
    id_ = 0;
  }
}

void
ShaderProgramARB::bind ()
{
  CHECK_OPENGL_ERROR();
  if (shaders_supported())
  {

    glLinkProgram(id_);

    // check to linking of the program
    GLint program_status[1];
    glGetProgramiv(id_, GL_LINK_STATUS, program_status);
    if (program_status[0] == GL_FALSE) {
      char program_log[1000];
      glGetInfoLogARB(id_, 1000, NULL, program_log);
      cerr << "Program Log: " << endl << program_log << endl;
    }

    glUseProgram(id_);

    /* get the variable and texture locations */
    const char *tex_strings[] = {"tex0", "tex1", "tex2", "tex3",
                                 "tex4", "tex5", "tex6", "tex7",
                                 "tex8", "tex9", "tex10", "tex11",
                                 "tex12", "tex13", "tex14", "tex15"};
    for (int i = 0; i < MAX_SHADER_UNIFORMS; i++) {
      int location = glGetUniformLocation(id_, tex_strings[i]);
      if (location != -1) { // able to get that link
  glUniform1i(location, i);
      }
    }
  }
  CHECK_OPENGL_ERROR();
}

void
ShaderProgramARB::release ()
{
  CHECK_OPENGL_ERROR();
  if (shaders_supported())
  {
    glUseProgram(0);
  }
  CHECK_OPENGL_ERROR();
}

int
ShaderProgramARB::activate()
{
  CHECK_OPENGL_ERROR();
  if (shaders_supported())
  {
    glUseProgram(id_);
  }
  CHECK_OPENGL_ERROR();

  return id_;
}


void
ShaderProgramARB::setTextures()
{
  CHECK_OPENGL_ERROR();
  if (shaders_supported())
  {
    for (int i = 0; i < MAX_SHADER_UNIFORMS; i++) {
      int location = glGetUniformLocation(id_, tex_strings[i]);
      if (location != -1) { // able to get that link
        glUniform1i(location, i);
      }
    }
  }
  CHECK_OPENGL_ERROR();
}

void
ShaderProgramARB::setLocalTexture(int i)
{
  CHECK_OPENGL_ERROR();
  if (shaders_supported())
  {
    int location = glGetUniformLocation(id_, tex_strings[i]);
    if (location != -1) // able to get that link
      glUniform1i(location, i);
  }
  CHECK_OPENGL_ERROR();
}

void
ShaderProgramARB::setLocalParam(int i, float x, float y, float z, float w)
{
  CHECK_OPENGL_ERROR();
  if (shaders_supported())
  {
    const char *loc_strings[] = {"loc0", "loc1", "loc2", "loc3",
                                 "loc4", "loc5", "loc6", "loc7",
                                 "loc8", "loc9", "loc10", "loc11",
                                 "loc12", "loc13", "loc14", "loc15"};
    
    int location = glGetUniformLocation(id_, loc_strings[i]);
    if (location == -1) {
      // cerr << "Error retrieving " << loc_strings[i] << " for " << i << endl;
      // cerr << program_ << endl;
    }
    else {
      glUniform4f(location, x, y, z, w);
    }
  }
  CHECK_OPENGL_ERROR();
}

VertexProgramARB::VertexProgramARB(const string& program)
  : ShaderProgramARB(program)
{
  type_ = GL_VERTEX_PROGRAM_ARB;
}

VertexProgramARB::~VertexProgramARB()
{
}

FragmentProgramARB::FragmentProgramARB(const string& program)
  : ShaderProgramARB(program)
{
  type_ = GL_FRAGMENT_PROGRAM_ARB;
}

FragmentProgramARB::~FragmentProgramARB()
{
}

} // end namespace SLIVR
