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
//    File   : TextureRenderer.cc
//    Author : Milan Ikits, Josh Stratton (changed to GLSL)
//    Date   : Wed Jul  7 23:34:20 2004


#include <slivr/gldefs.h>
#include <slivr/TextureRenderer.h>
#include <slivr/Color.h>

#include <slivr/Utils.h>
#include <slivr/ShaderProgramARB.h>
#include <slivr/VolShader.h>
#include <slivr/CM2Widget.h>
#include <slivr/CM2Shader.h>

#include <slivr/Colors.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using std::string;

namespace SLIVR {

//static bool mNV_float_buffer = true;

static const string Cmap2ShaderStringGLSL = 
"uniform vec4 loc0;\n"
"uniform sampler2D tex0;\n"
"\n"
"void main()\n"
"{\n"
"   vec4 c, z;\n"
"   vec4 t = gl_TexCoord[0];\n"
"   vec4 s = loc0; // {bp, sliceRatio, 0.0, 0.0}\n"
"   c = texture2D(tex0, t.st);\n"
"   z.w = pow(c.w, s.x); // pow(alpha, bp);\n"
"   z.w = 1.0 - z.w;\n"
"   c.w = pow(z.w, s.y);\n"
"   c.w = 1.0 - c.w;\n"
"   c.xyz = c.xyz * c.w;\n"
"   gl_FragColor = c;\n"
"}\n";

//#include <sstream>
//#include <string>
//inline std::string toStr(int x){ std::ostringstream ss; ss << x; std::string s(ss.str()); return s;}
void writeTextureToPPM( const char* fileName , GLuint tex, GLuint m_iSizeX, GLuint m_iSizeY, GLenum m_format, GLenum m_type)
{
    //cout << "m_iSizeX: " <<m_iSizeX << "  m_iSizeY: " << m_iSizeY << endl;
    GLuint m_iSizePerElement = sizeof(GL_UNSIGNED_BYTE);
    
    unsigned char*content = new unsigned char[m_iSizeX * m_iSizeY * m_iSizePerElement];
    glBindTexture(GL_TEXTURE_2D, tex);
    glGetTexImage(GL_TEXTURE_2D, 0, m_format, m_type, content);
    glBindTexture(GL_TEXTURE_2D, 0);
    
    Image img(m_iSizeX, m_iSizeY);
    // Image img(m_iSizeXX, m_iSizeYY);
    
    
    for(GLuint y = 0; y < m_iSizeY; y++){
        for(GLuint x = 0; x < m_iSizeX; x++){
            Colors C;
            size_t index = y * m_iSizeX + x;
            double alpha = double(content[4 * index + 3]) / 255.0;
            //alpha = 1.0;
            C.set_r(double(content[4 * index + 0]) / 255.0 * alpha);
            C.set_g(double(content[4 * index + 1]) / 255.0 * alpha);
            C.set_b(double(content[4 * index + 2]) / 255.0 * alpha);
            img.set(x, y, C);
        }
    }
    img.write(std::string(fileName) + ".ppm");
    
    delete [] content;
    content = NULL;
}


TextureRenderer::TextureRenderer(Texture* tex,
                                 ColorMap* cmap1, 
                                 const vector<ColorMap2*> &cmap2,
                                 int tex_mem) :
  tex_(tex),
  cmap1_(cmap1),
  cmap2_(cmap2),
  cmap1_dirty_(true),
  cmap2_dirty_(true),
  mode_(MODE_NONE),
  interp_(true),
  lighting_(false),
  sampling_rate_(1.0),
  irate_(0.5),
  imode_(false),
  slice_alpha_(0.5),
  cmap2_width_(256),
  cmap1_tex_(0),
  cmap2_tex_(0),
  tex_width_(0),
  blend_framebuffer_(0),
  blend_tex_id_(0),
  shader_factory_(0),
  cmap2_widget_framebuffer_(0),
  cmap2_widget_tex_id_(0),
  cmap2_framebuffer_(0),
  cmap2_tex_id_(0),
  cmap2_shader_glsl_(new FragmentProgramARB(Cmap2ShaderStringGLSL)),
  vol_shader_factory_(new VolShaderFactory()),
  blend_num_bits_(8),
  use_blend_buffer_(true),
  free_tex_mem_(tex_mem),
  use_stencil_(false),
  clear_pool_(false),
  reloadShaders_(true),
  readShadersFromFile_(true),
  textureDim_(256),
  shaderAlgo_(ALGO_NORMAL)
{
  attachmentpoints[0] = GL_COLOR_ATTACHMENT0_EXT;
  attachmentpoints[1] = GL_COLOR_ATTACHMENT1_EXT;
  attachmentpoints[2] = GL_COLOR_ATTACHMENT2_EXT;
  attachmentpoints[3] = GL_COLOR_ATTACHMENT3_EXT;
  
  attachmentpoints1[0] = GL_COLOR_ATTACHMENT0_EXT;
  attachmentpoints1[1] = GL_COLOR_ATTACHMENT2_EXT;
  attachmentpoints2[0] = GL_COLOR_ATTACHMENT1_EXT;
  attachmentpoints2[1] = GL_COLOR_ATTACHMENT3_EXT;
}

TextureRenderer::TextureRenderer(const TextureRenderer& copy) :
  tex_(copy.tex_),
  cmap1_(copy.cmap1_),
  cmap2_(copy.cmap2_),
  cmap1_dirty_(copy.cmap1_dirty_),
  cmap2_dirty_(copy.cmap2_dirty_),
  mode_(copy.mode_),
  interp_(copy.interp_),
  lighting_(copy.lighting_),
  sampling_rate_(copy.sampling_rate_),
  irate_(copy.irate_),
  imode_(copy.imode_),
  slice_alpha_(copy.slice_alpha_),
  cmap2_width_(copy.cmap2_width_),
  cmap1_tex_(copy.cmap1_tex_),
  cmap2_tex_(copy.cmap2_tex_),
  shader_factory_(copy.shader_factory_),
  cmap2_shader_glsl_(copy.cmap2_shader_glsl_),
  vol_shader_factory_(copy.vol_shader_factory_),
  blend_num_bits_(copy.blend_num_bits_),
  use_blend_buffer_(copy.use_blend_buffer_),
  free_tex_mem_(copy.free_tex_mem_),
  use_stencil_(copy.use_stencil_),
  clear_pool_(true),
  reloadShaders_(true),
  readShadersFromFile_(true),
  textureDim_(256),
  shaderAlgo_(ALGO_NORMAL)
{
  attachmentpoints[0] = GL_COLOR_ATTACHMENT0_EXT;
  attachmentpoints[1] = GL_COLOR_ATTACHMENT1_EXT;
  attachmentpoints[2] = GL_COLOR_ATTACHMENT2_EXT;
  attachmentpoints[3] = GL_COLOR_ATTACHMENT3_EXT;
  
  attachmentpoints1[0] = GL_COLOR_ATTACHMENT0_EXT;
  attachmentpoints1[1] = GL_COLOR_ATTACHMENT2_EXT;
  attachmentpoints2[0] = GL_COLOR_ATTACHMENT1_EXT;
  attachmentpoints2[1] = GL_COLOR_ATTACHMENT3_EXT;
}


TextureRenderer::~TextureRenderer()
{
  delete cmap2_shader_glsl_;
  delete vol_shader_factory_;
  
  for (size_t i = 0; i < tex_pool_.size(); i++)
  {
    glDeleteTextures(1, (GLuint*)&tex_pool_[i].id);  
  }

  glDeleteFramebuffersEXT(1, &blend_framebuffer_);
  glDeleteTextures(1, &blend_tex_id_);
}


void
TextureRenderer::set_texture(Texture* tex)
{
  if (tex_ != tex) 
  {
    tex_ = tex;
    // new texture, flag existing tex id's for deletion.
    clear_pool_ = true;
  }
}

void
TextureRenderer::set_colormap1(ColorMap* cmap1)
{
  cmap1_ = cmap1;
  cmap1_dirty_ = true;
}

void
TextureRenderer::set_colormap2(const vector<ColorMap2*> &cmap2)
{
  cmap2_ = cmap2;
  cmap2_dirty_ = true;
}

void
TextureRenderer::set_colormap2_width(int size)
{
  if (cmap2_width_ != size) {
    cmap2_width_ = size;
    cmap2_dirty_ = true;
  }
}

void
TextureRenderer::set_slice_alpha(double alpha)
{
  if(fabs(slice_alpha_ - alpha) > 0.0001) {
    slice_alpha_ = alpha;
    alpha_dirty_ = true;
  }
}


void
TextureRenderer::set_blend_num_bits(int b)
{
  blend_num_bits_ = b;
}

bool
TextureRenderer::use_blend_buffer()
{
  return use_blend_buffer_;
}

Ray
TextureRenderer::compute_view()
{
  const Transform &field_trans = tex_->transform();
  double mvmat[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, mvmat);
  // index space view direction
  Vector v = field_trans.unproject(Vector(-mvmat[2], -mvmat[6], -mvmat[10]));
  v.safe_normalize();
  Transform mv;
  mv.set_trans(mvmat);
  Point p = field_trans.unproject(mv.unproject(Point(0,0,0)));
  return Ray(p, v);
}

bool
TextureRenderer::test_against_view(const BBox &bbox)
{
  double mvmat[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, mvmat);
  Transform mv;
  mv.set_trans(mvmat);

  double prmat[16];
  glGetDoublev(GL_PROJECTION_MATRIX, prmat);
  Transform pr;
  pr.set_trans(prmat);

  bool overx = true;
  bool overy = true;
  bool overz = true;
  bool underx = true;
  bool undery = true;
  bool underz = true;
  for (int i = 0; i < 8; i++)
  {
    const Point pold((i&1)?bbox.min().x():bbox.max().x(),
                     (i&2)?bbox.min().y():bbox.max().y(),
                     (i&4)?bbox.min().z():bbox.max().z());
    const Point p = pr.project(mv.project(pold));
    overx = overx && (p.x() > 1.0);
    overy = overy && (p.y() > 1.0);
    overz = overz && (p.z() > 1.0);
    underx = underx && (p.x() < -1.0);
    undery = undery && (p.y() < -1.0);
    underz = underz && (p.z() < -1.0);
  }

  return !(overx || overy || overz || underx || undery || underz);
}


// Pool applies only to the current texture. 
// Clear it with each new set_texture.
void
TextureRenderer::clear_tex_pool() 
{
  for(unsigned int i = 0; i < tex_pool_.size(); i++)
  {
    size_t size = (tex_pool_[i].nx * tex_pool_[i].ny * 
       tex_pool_[i].nz * tex_pool_[i].nb);
    // delete tex object.
    glDeleteTextures(1, (GLuint*)&tex_pool_[i].id);
    tex_pool_[i].id = 0;
    free_tex_mem_ += size;
  }
  clear_pool_ = false;
}


void
TextureRenderer::load_brick(vector<TextureBrick*> &bricks, int bindex,
                            bool use_cmap2)
{
  CHECK_OPENGL_ERROR();
  if (clear_pool_) clear_tex_pool();
  TextureBrick* brick = bricks[bindex];
  int nc = use_cmap2 ? brick->nc() : 1;
#if !defined(GL_ARB_fragment_program) && !defined(GL_ATI_fragment_shader)
  nc = 1;
#endif
  int idx[2];
  for(int c = 0; c < nc; c++) 
  {
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
    if (glActiveTexture) { glActiveTexture(GL_TEXTURE0+c); }
#endif
    int nb = brick->nb(c);
    int nx = brick->nx();
    int ny = brick->ny();
    int nz = brick->nz();
    GLenum textype = brick->tex_type(c);

    //! Try to find the existing texture in tex_pool_, for this brick.
    idx[c] = -1;
    for(unsigned int i = 0; i < tex_pool_.size() && idx[c] < 0; i++)
    {
      if(tex_pool_[i].id != 0 && tex_pool_[i].brick == brick
         && !brick->dirty() && tex_pool_[i].comp == c
         && nx == tex_pool_[i].nx && ny == tex_pool_[i].ny
         && nz == tex_pool_[i].nz && nb == tex_pool_[i].nb
         && textype == tex_pool_[i].textype
         && glIsTexture(tex_pool_[i].id))
      {
        if (tex_pool_[i].brick == brick)
        {
          idx[c] = i;
        }
        else
        {
          bool found = false;
          for (unsigned int j = 0; j < bricks.size(); j++)
          {
            if (bricks[j] == brick)
            {
              found = true;
            }
          }
          if (!found)
          {
            idx[c] = i;
          }
        }
      }
    }
    if(idx[c] != -1) 
    {
      //! The texture object was located, bind it.
      // bind texture object
      glBindTexture(GL_TEXTURE_3D, tex_pool_[idx[c]].id);
      // set interpolation method
      if(interp_) {
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      } else {
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      }
    } else {
      //! try again to find the matching texture object
      for(unsigned int i=0; i<tex_pool_.size() && idx[c]<0; i++)
      {
        if(tex_pool_[i].id != 0 && c == tex_pool_[i].comp
           && nx == tex_pool_[i].nx && ny == tex_pool_[i].ny
           && nz == tex_pool_[i].nz && nb == tex_pool_[i].nb
           && textype == tex_pool_[i].textype)
        {
          if (tex_pool_[i].brick == brick)
          {
            idx[c] = i;
          }
          else
          {
            bool found = false;
            for (unsigned int j = 0; j < bricks.size(); j++)
            {
              if (bricks[j] == brick)
              {
                found = true;
              }
            }
            if (!found)
            {
              idx[c] = i;
            }
          }
        }
      }
      bool reuse = (idx[c] >= 0);
      // 
      if(!reuse) {
        // if cannot reuse existing object allocate new object
        int new_size = nx*ny*nz*nb;
        if(new_size > free_tex_mem_) {
          // if there's no space, find object to replace
          int free_idx = -1;
          int size, size_max = -1;
          // find smallest available objects to delete
          // TODO: this is pretty dumb, optimize it later
          for(unsigned int i=0; i<tex_pool_.size(); i++) {
            for(int j=0; j<c; j++) {
              if(idx[j] == (int)i) continue;
            }
            size = tex_pool_[i].nx*tex_pool_[i].ny*tex_pool_[i].nz*tex_pool_[i].nb;
            if(new_size < free_tex_mem_+size && (size_max < 0 || size < size_max)) {
              free_idx = i;
              size_max = size;
            }
          }
          if (free_idx != -1)
          {
            // delete found object            
            glDeleteTextures(1, (GLuint*)&tex_pool_[free_idx].id);
            tex_pool_[free_idx].id = 0;
          }
          free_tex_mem_ += size_max;
        }
        // find tex table entry to reuse
        for(unsigned int i=0; i<tex_pool_.size() && idx[c]<0; i++) {
          if(tex_pool_[i].id == 0)
            idx[c] = i;
        }
        // allocate new object
        unsigned int tex_id;
        glGenTextures(1, (GLuint*)&tex_id);
        if(idx[c] < 0) {
          // create new entry
          tex_pool_.push_back(TexParam(nx, ny, nz, nb, textype, tex_id));
          idx[c] = tex_pool_.size()-1;
        } else {
          // reuse existing entry
          tex_pool_[idx[c]].nx = nx; tex_pool_[idx[c]].ny = ny;
          tex_pool_[idx[c]].nz = nz; tex_pool_[idx[c]].nb = nb;
          //make sure the previous texture has been deleted
          glDeleteTextures(1, (GLuint*)&tex_pool_[idx[c]].id);
          tex_pool_[idx[c]].id = tex_id;
        }
        free_tex_mem_ -= new_size;
      }
      tex_pool_[idx[c]].brick = brick;
      tex_pool_[idx[c]].comp = c;
      // bind texture object
      glPixelStorei(GL_UNPACK_ALIGNMENT,nb);
      glBindTexture(GL_TEXTURE_3D, tex_pool_[idx[c]].id);
      // set border behavior
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
      // set interpolation method
      if(interp_) {
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      } else {
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      }

      // download texture data
      glPixelStorei(GL_UNPACK_ROW_LENGTH, brick->sx());
      glPixelStorei(GL_UNPACK_IMAGE_HEIGHT, brick->sy());
      glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);    // needed on osx mtn lion
      glPixelStorei(GL_UNPACK_ALIGNMENT, (nb == 1)?1:4);
#if defined( GL_TEXTURE_COLOR_TABLE_SGI ) && defined(__sgi)
      if (reuse)
      {
        glTexSubImage3DEXT(GL_TEXTURE_3D, 0, 0, 0, 0, nx, ny, nz, GL_RED,
                           brick->tex_type(c), brick->tex_data(c));
      }
      else
      {
        glTexImage3DEXT(GL_TEXTURE_3D, 0, GL_INTENSITY8,
                        nx, ny, nz, 0, GL_RED,
                        brick->tex_type(c), brick->tex_data(c));
      }
#elif defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
      if (ShaderProgramARB::shaders_supported())
      {
        unsigned int format = (nb == 1 ? GL_LUMINANCE : GL_RGBA);
        if (reuse && glTexSubImage3D)
        {
          glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, nx, ny, nz, format,
                          brick->tex_type(c), brick->tex_data(c));
        }
        else
        {
          if (glTexImage3D)
            glTexImage3D(GL_TEXTURE_3D, 0, format, nx, ny, nz, 0, format,
                         brick->tex_type(c), brick->tex_data(c));
        }
      }
      else
#endif
#if !defined(__sgi)
#  if defined(GL_EXT_shared_texture_palette) && !defined(__APPLE__)
      {
        if (reuse && glTexSubImage3D)
        {
          glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, nx, ny, nz,
                          GL_COLOR_INDEX,
                          brick->tex_type(c), brick->tex_data(c));
        }
        else
        {
          if (glTexImage3D)
            glTexImage3D(GL_TEXTURE_3D, 0, GL_COLOR_INDEX8_EXT,
                         nx, ny, nz, 0, GL_COLOR_INDEX,
                         brick->tex_type(c), brick->tex_data(c));
        }
      }
#  elif defined(_WIN32) || defined(GL_VERSION_1_2) // Workaround for old bad nvidia headers.
      {
        if (reuse && glTexSubImage3D)
        {
          glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, nx, ny, nz,
                          GL_LUMINANCE,
                          brick->tex_type(c), brick->tex_data(c));
        }
        else
        {
          if (glTexImage3D) {
            glTexImage3D(GL_TEXTURE_3D, 0, GL_LUMINANCE,
                         nx, ny, nz, 0, GL_LUMINANCE,
                         brick->tex_type(c), brick->tex_data(c));
    }
  }
      }
#  endif
#endif // !__sgi
      glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
      glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
      glPixelStorei(GL_UNPACK_IMAGE_HEIGHT, 0);
      glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    }
  }
  brick->set_dirty(false);
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
  if (glActiveTexture) { glActiveTexture(GL_TEXTURE0); }
#endif
  CHECK_OPENGL_ERROR();
}

void
TextureRenderer::draw_polygons(vector<float>& vertex, 
                               vector<float>& texcoord,
                               vector<int>& poly, 
                               bool normal, bool fog,
                               vector<int> *mask, 
                               ShaderProgramARB* shader)

{
  CHECK_OPENGL_ERROR();

  float mvmat[16];
  if(fog) {
    glGetFloatv(GL_MODELVIEW_MATRIX, mvmat);
  }

  for(unsigned int i=0, k=0; i<poly.size(); i++) {
    if (mask && shader) {
      int maskval = ((*mask)[i] << 1) + 1;
      if (!maskval) {
        k += poly[i];
        continue;
      }
      const float v = float(maskval);
      shader->setLocalParam(3, v, v, v, v);
    }

    glBegin(GL_POLYGON);
    {
      if(normal) {
        float* v0 = &vertex[(k+0)*3];
        float* v1 = &vertex[(k+1)*3];
        float* v2 = &vertex[(k+2)*3];
        Vector dv1(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
        Vector dv2(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);
        Vector n = Cross(dv1, dv2);
        n.normalize();
        glNormal3f(n.x(), n.y(), n.z());
      }
      for(int j=0; j<poly[i]; j++) {
        float* t = &texcoord[(k+j)*3];
        float* v = &vertex[(k+j)*3];
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
#  ifdef _WIN32
        if (glMultiTexCoord3f) {
#  endif // _WIN32
          glMultiTexCoord3f(GL_TEXTURE0, t[0], t[1], t[2]);
          if(fog) {
            float vz = mvmat[2]*v[0] + mvmat[6]*v[1] + mvmat[10]*v[2] + mvmat[14];
            glMultiTexCoord3f(GL_TEXTURE1, -vz, 0.0, 0.0);
          }
#  ifdef _WIN32
        } else {
          glTexCoord3f(t[0], t[1], t[2]);
        }
#  endif // _WIN32
#else
        glTexCoord3f(t[0], t[1], t[2]);
#endif
        glVertex3f(v[0], v[1], v[2]);
      }
    }
    glEnd();

    k += poly[i];
  }
  CHECK_OPENGL_ERROR();
}

void
TextureRenderer::draw_polygons_wireframe(vector<float>& vertex,
                                         vector<float>& /*texcoord*/,
                                         vector<int>& poly,
                                         bool /*normal*/, bool fog,
                                         vector<int> *mask)
{
  //  di_->polycount_ += poly.size();
  float mvmat[16];
  if(fog) {
    glGetFloatv(GL_MODELVIEW_MATRIX, mvmat);
  }
  for(unsigned int i=0, k=0; i<poly.size(); i++)
  {
    if (mask) {
      int v = (*mask)[i] ? 2:1;
      glColor4d(v & 1 ? 1.0 : 0.0, 
                v & 2 ? 1.0 : 0.0, 
                v & 4 ? 1.0 : 0.0, 1.0);
    }
      
    glBegin(GL_LINE_LOOP);
    {
      for(int j=0; j<poly[i]; j++)
      {
        float* v = &vertex[(k+j)*3];
        glVertex3f(v[0], v[1], v[2]);
      }
    }
    glEnd();
    k += poly[i];
  }
}


void
TextureRenderer::build_colormap1(vector<float> cmap_array,
                                 unsigned int& cmap_tex, bool& cmap_dirty,
                                 bool& alpha_dirty,  double level_exponent)
{
  if(cmap_dirty || alpha_dirty) {
    bool size_dirty = false;
    if(1024 != cmap_array.size()) {
      cmap_array.resize(1024);
      size_dirty = true;
    }
    // rebuild texture
    double dv = 1.0 / 256; // colors defined from (c >= 0.0, c < 1.0)
    switch(mode_) {
    case MODE_SLICE: {
      for(int j=0; j < 1024; j+=4) {
        // interpolate from colormap
  double t = (j/4)*dv;
  float r,g,b,alpha;
  cmap1_->get_color(t, r, g, b, alpha);
        // pre-multiply and quantize
        cmap_array[j+0] = r*alpha;
        cmap_array[j+1] = g*alpha;
        cmap_array[j+2] = b*alpha;
        cmap_array[j+3] = alpha;
      }
    } break;
    case MODE_MIP: {
      for(int j=0; j < 1024; j+=4) {
        // interpolate from colormap
  double t = (j/4)*dv;
  float r,g,b,alpha;
  cmap1_->get_color(t, r, g, b, alpha);
        // pre-multiply and quantize
        cmap_array[j+0] = r*alpha;
        cmap_array[j+1] = g*alpha;
        cmap_array[j+2] = b*alpha;
        cmap_array[j+3] = alpha;
      }
    } break;
    case MODE_OVER: {
      double bp = tan(1.570796327 * (0.5 - slice_alpha_*0.49999));
      for(int j=0; j < 1024; j+=4) {
        // interpolate from colormap
        double t = (j/4)*dv;
        float r,g,b,alpha;
        cmap1_->get_color(t, r, g, b, alpha);
        // scale slice opacity
        alpha = pow(alpha, (float)bp);
        // opacity correction
        //alpha = 1.0 - pow((1.0 - alpha), 1.0/sampling_rate);
        alpha = 1.0 - pow((1.0 - alpha), imode_ ?
                          1.0/irate_/pow(2.0, level_exponent) :
                          1.0/sampling_rate_/pow(2.0, level_exponent) );

        // pre-multiply and quantize
        cmap_array[j+0] = r*alpha;
        cmap_array[j+1] = g*alpha;
        cmap_array[j+2] = b*alpha;
        cmap_array[j+3] = alpha;
      }
    } break;
    default:
      break;
    }
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
    // This texture is not used if there is no shaders.
    // glColorTable is used instead.
    if (ShaderProgramARB::shaders_supported())
    {
      // Update 1D texture.
      if (cmap_tex == 0 || size_dirty)
      {
        glDeleteTextures(1, (GLuint*)&cmap_tex);
        glGenTextures(1, (GLuint*)&cmap_tex);
        glBindTexture(GL_TEXTURE_1D, cmap_tex);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA16, 256, 0,
                     GL_RGBA, GL_FLOAT, &cmap_array[0]);
      }
      else
      {
        glBindTexture(GL_TEXTURE_1D, cmap_tex);
        glTexSubImage1D(GL_TEXTURE_1D, 0, 0, 256,
                        GL_RGBA, GL_FLOAT, &cmap_array[0]);
      }
    }
#endif
    cmap_dirty = false;
    alpha_dirty = false;
  }
}


void
TextureRenderer::colormap2_hardware_rasterize() 
{
  CHECK_OPENGL_ERROR();
  if (cmap2_dirty_ || alpha_dirty_)
  {
    // Cache this value to reset, in case another framebuffer is active.
    GLint cur_framebuffer_id;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING_EXT, &cur_framebuffer_id);
    GLint cur_draw_buffer;
    glGetIntegerv(GL_DRAW_BUFFER, &cur_draw_buffer);
    GLint cur_read_buffer;
    glGetIntegerv(GL_READ_BUFFER, &cur_read_buffer);

    glPushAttrib(GL_VIEWPORT_BIT | GL_ENABLE_BIT);
    glViewport(0, 0, tex_width_, 64);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(-1.0, -1.0, 0.0);
    glScalef(2.0/Pow2(cmap2_.size()), 2.0, 2.0);

    glDisable(GL_DEPTH_TEST);  // GL_ENABLE_BIT
    glDisable(GL_LIGHTING);    // GL_ENABLE_BIT
    glDisable(GL_CULL_FACE);   // GL_ENABLE_BIT

    GLfloat clear_color[4];
    glGetFloatv(GL_COLOR_CLEAR_VALUE, clear_color);

    // Disable clipping planes.  Reset with GL_ENABLE_BIT attrib
    for (int num = 0; num < 6; num++) {
      glDisable((GLenum)(GL_CLIP_PLANE0+num));
    }

    //--------------------------------------------------------------
    // Do the widget draw
    if (cmap2_dirty_)
    {
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
      if (glActiveTexture) { glActiveTexture(GL_TEXTURE0); }
#endif

      // Set up the widget framebuffer
      if (cmap2_widget_framebuffer_ == 0)
      {
        glGenFramebuffersEXT(1, &cmap2_widget_framebuffer_);
      }
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, cmap2_widget_framebuffer_);
    
      // Set up the widget texture.
      if (cmap2_widget_tex_id_ == 0)
      {
        glGenTextures(1, &cmap2_widget_tex_id_);
        glBindTexture(GL_TEXTURE_2D, cmap2_widget_tex_id_);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tex_width_, 64, 0,
                     GL_RGBA, GL_INT, NULL);
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                                  GL_COLOR_ATTACHMENT0_EXT,
                                  GL_TEXTURE_2D, cmap2_widget_tex_id_, 0);

        CHECK_FRAMEBUFFER_STATUS();
      }

      glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
      glClearColor(0.0, 0.0, 0.0, 0.0);
      glClear(GL_COLOR_BUFFER_BIT);

      glEnable(GL_BLEND);
      // blend mode?

      // Rasterize widgets
      for (unsigned int c = 0; c < cmap2_.size(); ++c)
      {
        vector<CM2Widget*> &widgets = cmap2_[c]->widgets();
        for (unsigned int i=0; i<widgets.size(); i++)
        {
          // Render to the texture memory.
          widgets[i]->rasterize(*shader_factory_);
          CHECK_OPENGL_ERROR();
        }
        // Translate to the next colormap location.
        glTranslatef(1.0, 0.0, 0.0);
      }
      alpha_dirty_ = true;
    }

    //--------------------------------------------------------------
    // Do the alpha draw

    // Unshift the multiple CM2 maps for the alpha blend.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-1.0, -1.0, 0.0);
    glScalef(2.0, 2.0, 2.0);

    //    GLint val;
    //    glGetIntegerv( GL_MAX_COLOR_ATTACHMENTS_EXT, &val );
    //    cout << "MAX_COLOR_ATTACHMENTS_EXT is " << val << "\n";

#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
    if (glActiveTexture) { glActiveTexture(GL_TEXTURE0); }
#endif

    // Set up the framebuffer
    if (cmap2_framebuffer_ == 0)
    {
      glGenFramebuffersEXT(1, &cmap2_framebuffer_);
    }
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, cmap2_framebuffer_);

    // Set up the colormap2 texture.
    if (cmap2_tex_id_ == 0)
    {
      glGenTextures(1, &cmap2_tex_id_);
      glBindTexture(GL_TEXTURE_2D, cmap2_tex_id_);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tex_width_, 64, 0,
                   GL_RGBA, GL_INT, NULL);

      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
        GL_COLOR_ATTACHMENT0_EXT,
        GL_TEXTURE_2D, cmap2_tex_id_, 0);

      CHECK_FRAMEBUFFER_STATUS();
    }

    glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glDisable(GL_BLEND);

    // Draw the widget texture with the alpha computing shader.
    FragmentProgramARB* shader = cmap2_shader_glsl_;
    shader->bind();
    double bp = mode_ == MODE_MIP ? 1.0 : 
      tan(1.570796327 * (0.5 - slice_alpha_*0.49999));
    shader->setLocalParam(0, bp, imode_ ? 1.0/irate_ : 1.0/sampling_rate_, 
                          0.0, 0.0);
    CHECK_OPENGL_ERROR();

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, cmap2_widget_tex_id_);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glBegin(GL_QUADS);
    {
      glTexCoord2f( 0.0,  0.0);
      glVertex2f( 0.0,  0.0);
      glTexCoord2f(1.0,  0.0);
      glVertex2f( 1.0,  0.0);
      glTexCoord2f(1.0,  1.0);
      glVertex2f( 1.0,  1.0);
      glTexCoord2f( 0.0,  1.0);
      glVertex2f( 0.0,  1.0);
    }
    glEnd();
    shader->release();
    CHECK_OPENGL_ERROR();    

    // Clean up the OpenGL state.
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glPopAttrib();
    glClearColor(clear_color[0], clear_color[1], 
     clear_color[2], clear_color[3]);

    // Restore buffer state.
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, cur_framebuffer_id); 
    glDrawBuffer(cur_draw_buffer);
    glReadBuffer(cur_read_buffer);
  }
}


void 
TextureRenderer::initEyeTex(int texWidth, int texHeight, int value1=0, int value2=255){
  unsigned char *textureArray = new unsigned char[4*texWidth*texHeight];
  for (int i=0; i<(textureDim_*textureDim_); i++)
    textureArray[i*4 +0] = textureArray[i*4 +1] = textureArray[i*4 +2] =  textureArray[i*4 +3] = value1;

    // Eye buffer
    glEnable(GL_TEXTURE_2D);
    for (int i=0; i<2; i++){
        glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id0[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, texWidth, texHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, textureArray);
    }
    
    // Occlusion Buffer - initialize to white
    for (int i=0; i<(textureDim_*textureDim_); i++)
    textureArray[i*4 +0] = textureArray[i*4 +1] = textureArray[i*4 +2] =  textureArray[i*4 +3] = value2;
    
    for (int i=0; i<2; i++){
        glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id1[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, texWidth, texHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, textureArray);
    }
    
    delete []textureArray;
}


void
TextureRenderer::init_textures()
{
  glGenTextures(2, eyeBuffer_tex_id0);
    for (int i=0; i<2; i++){
        glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id0[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, textureDim_, textureDim_, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    }

  glGenTextures(2, eyeBuffer_tex_id1);
    for (int i=0; i<2; i++){
        glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id1[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, textureDim_, textureDim_, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    }
}

void
TextureRenderer::destroy_textures()
{
  glDeleteTextures(2, eyeBuffer_tex_id0);
  glDeleteTextures(2, eyeBuffer_tex_id1);
}


bool
TextureRenderer::init_FrameBuffer()
{
  glGenFramebuffersEXT(1, &frmBuffer_);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frmBuffer_);

  // attach two textures to FBO
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachmentpoints[0], GL_TEXTURE_2D, eyeBuffer_tex_id0[0], 0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachmentpoints[1], GL_TEXTURE_2D, eyeBuffer_tex_id0[1], 0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachmentpoints[2], GL_TEXTURE_2D, eyeBuffer_tex_id1[0], 0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachmentpoints[3], GL_TEXTURE_2D, eyeBuffer_tex_id1[1], 0);

  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  if (status == GL_FRAMEBUFFER_COMPLETE_EXT){
  //  cout << "All Good!!!" << endl;
  }
  else{
    cout << "Error: Binding framebuffer failed" << status << endl;
    return false;
  }
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  return true;
}


void
TextureRenderer::destroy_FrameBuffer()
{
  glDeleteFramebuffers(1, &frmBuffer_);
}


void 
TextureRenderer::draw_polygonsOccSh(vector<float>& vertex,
                               vector<float>& texcoord,
                               vector<int>& poly,
                               bool normal, bool fog,
                               vector<int> *mask,
                               ShaderProgramARB *shaderOccSh, ShaderProgramARB *shaderTexDisp)
{
  CHECK_OPENGL_ERROR();
  float mvmat[16];
  if (fog) {
    glGetFloatv(GL_MODELVIEW_MATRIX, mvmat);
  }

  //int programId, location;
  int location;
  int writeTex, readTex;
  int lastWritePhase1, lastWritePhase2;
  lastWritePhase1 = lastWritePhase2 = 0;

  int z;          // current z slice position
  int totalVerticesInPoly = 0;
  for (unsigned int i=0; i<poly.size(); i++)
    totalVerticesInPoly += poly[i];
    
  float lengthSlice = textureDim_/(float)poly.size();
  float lengthOpp = (lengthSlice * tan(occSh_blurAngle_ * PI/180.0))/2.0;
    
  writeTex = 1; readTex = 0;

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glPushMatrix();

  glDisable(GL_BLEND);
    initEyeTex(textureDim_,textureDim_);

    // First Phase: Render to eyebuffer
    glViewport(0,0,textureDim_,textureDim_);  // set the viewport for now to the texture size
    screenSize[0] = screenSize[1] = textureDim_;

    int programId = shaderOccSh->activate();
    location = glGetUniformLocation(programId, "windowSize");
    glUniform2fv(location, 1, screenSize);

    location = glGetUniformLocation(programId, "angleFactor");
    glUniform1fARB(location, lengthOpp);

    location = glGetUniformLocation(programId, "ambIntensity");
    glUniform1fARB(location, occSh_ambIntensity_);

 
    z = -1;
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frmBuffer_);
    for (unsigned int i=poly.size()-1, k=totalVerticesInPoly; i>=1; i--) {  

      if (writeTex == 0)  // need to write to eye buffer + to occlusion buffer
          glDrawBuffers(2, attachmentpoints1);
      else
          glDrawBuffers(2, attachmentpoints2);
      z++;
      
      glEnable(GL_TEXTURE_2D);
      // Eye Buffer
      glActiveTexture(GL_TEXTURE6);
      glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id0[readTex]);
      shaderOccSh->setLocalTexture(6);
      
      // Occlusion Buffer
      glActiveTexture(GL_TEXTURE7);
      glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id1[readTex]);
      shaderOccSh->setLocalTexture(7);

      glBegin(GL_POLYGON);
      {
        if (normal) {
          float* v0 = &vertex[(k+0)*3];
          float* v1 = &vertex[(k+1)*3];
          float* v2 = &vertex[(k+2)*3];
          Vector dv1(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
          Vector dv2(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);
          Vector n = Cross(dv1, dv2);   n.normalize();
          glNormal3f(n.x(), n.y(), n.z());
        }

        // draw the number of vertices making up the polygon
        for(int j=0; j<poly[i]; j++) {
          float* t = &texcoord[(k+j)*3];
          float* v = &vertex[(k+j)*3];

          #if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
          #ifdef _WIN32
          if (glMultiTexCoord3f) {
          #endif // _WIN32
            glMultiTexCoord3f(GL_TEXTURE0, t[0],   t[1],     t[2]);
            if (fog) {
              float vz = mvmat[2]*v[0] + mvmat[6]*v[1] + mvmat[10]*v[2] + mvmat[14];
              glMultiTexCoord3f(GL_TEXTURE1, -vz, 0.0, 0.0);
            }
          #ifdef _WIN32
          }
          else{
            glTexCoord3f(t[0], t[1], t[2]);
          }
          #endif // _WIN32
          #else
          glTexCoord3f(t[0], t[1], t[2]);
          #endif
          glVertex3f(v[0], v[1], v[2]);
        }
      }
      glEnd();

      // swap read and write textures
      writeTex = 1 - writeTex;
      readTex = 1 - readTex;

      k -= poly[i]; 
    }
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      
    lastWritePhase1 = readTex;
    shaderOccSh->release();

    glEnable(GL_BLEND);


    // Final Phase: Render to texture
    glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
    glPushMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
        glLoadIdentity();
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
          glLoadIdentity();
          programId = shaderTexDisp->activate();
            glEnable(GL_TEXTURE_2D);

            location = glGetUniformLocation(programId, "tex7");
            glActiveTexture(GL_TEXTURE7);
            glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id0[lastWritePhase1]);
            glUniform1i(location, 7);
            
            location = glGetUniformLocation(programId, "option");
            glUniform1i(location, 0);

            glBegin (GL_QUADS);
              glTexCoord2f(0.0, 0.0);  glVertex3f (-1, -1, 0.0);
              glTexCoord2f(1.0, 0.0);  glVertex3f ( 1, -1, 0.0);
              glTexCoord2f(1.0, 1.0);  glVertex3f ( 1,  1, 0.0);
              glTexCoord2f(0.0, 1.0);  glVertex3f (-1,  1, 0.0);
            glEnd ();

          shaderTexDisp->release();
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
    glPopMatrix();

  glPopMatrix();
  glPopAttrib();

  CHECK_OPENGL_ERROR();
}





void TextureRenderer::draw_polygonsDOF(vector<float>& vertex,
                               vector<float>& texcoord,
                               vector<int>& poly,
                               bool normal, bool fog,
                               vector<int> *mask,
                               ShaderProgramARB *shaderDoF, ShaderProgramARB *shaderTexDisp)
{
  CHECK_OPENGL_ERROR();
  float mvmat[16];
  if (fog) {
    glGetFloatv(GL_MODELVIEW_MATRIX, mvmat);
  }

  int location;
  int writeTex, readTex;
  int lastWritePhase1, lastWritePhase2;
  lastWritePhase1 = lastWritePhase2 = 0;

  float Cz, radS;
  int z;          // current z slice position
  int zf;         // focal plane
  int zf_slicer;  // focal plane for slicer operation order - from poly.size to 0 for front to back

  zf = poly.size() * dof_focusPosition_;
  zf_slicer = poly.size() - (poly.size() * dof_focusPosition_);
      
  int totalVerticesInPoly = 0;
  for (unsigned int i=0; i<poly.size(); i++)
    totalVerticesInPoly += poly[i];
    
  float focusRange = (dof_focusRange_ * poly.size());
  float lengthSlice = textureDim_/(float)poly.size();
  float lengthOpp = (lengthSlice * tan(dof_blurAngle_ * PI/180.0))/2.0;
  screenSize[0] = screenSize[1] = textureDim_;


  cout << "-------------------------------------------" << endl;
  cout << "dof_focusPosition_: " << dof_focusPosition_ << endl;
  cout << "zf: " << zf << endl;
  cout << "zf_slicer: " << zf_slicer << endl;
  cout << "poly.size(): " << poly.size() << endl;
  cout << "focusRange: " << focusRange << endl;
  cout << "focusRange: " << focusRange << endl;
  cout << "lengthSlice: " << lengthSlice << endl;
  cout << "lengthOpp: " << lengthOpp << endl;
  cout << "dof_focusMode_: " << dof_focusMode_ << endl;
  cout << "screenSize/textureDim_: " << screenSize[0] << ", " << screenSize[1] << endl;
  if (dof_focusMode_ == DOF_USER)
    cout << "DOF USER" << endl;
  else
    cout << "DOF AUTO" << endl;

  writeTex = 1; readTex = 0;

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glPushMatrix();

    glDisable(GL_BLEND);
    initEyeTex(textureDim_,textureDim_,0,0);

    glViewport(0,0,textureDim_,textureDim_);
    
    int programId = shaderDoF->activate();

    location = glGetUniformLocation(programId, "windowSize");
    glUniform2fv(location, 1, screenSize);

    location = glGetUniformLocation(programId, "angleFactor");
    glUniform1fARB(location, lengthOpp);

             
  // Phase 1: front to in-focus-plane (front-to-back)
    location = glGetUniformLocation(programId, "traversalDirection");
    glUniform1i(location, FRONT_TO_BACK);

    z = 0;
    for (unsigned int i=poly.size()-1, k=totalVerticesInPoly; i>=zf_slicer; i--){  
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frmBuffer_);
      glDrawBuffer(attachmentpoints[writeTex]);
      
      z++;
      if (dof_focusMode_ == DOF_USER){
        if ( (z >= (zf-focusRange)) && (z <= (zf+focusRange)) )
            radS = 1.0;
        else
            radS = 4.0;
      }
      else{
        Cz = dof_aperture_ * ((zf-z)/(float)z);

        if (radS <= dof_threshold_)
            radS = 1.0;
        else
            radS = Cz/dof_threshold_;
      }

              
      if (radS <= 2.0){
        location = glGetUniformLocation(programId, "numSamples");
        glUniform1i(location, 1);
      }
      else
        if (radS <= 5.0){
          location = glGetUniformLocation(programId, "numSamples");
          glUniform1i(location, 5);
        }
        else{
          location = glGetUniformLocation(programId, "numSamples");
          glUniform1i(location, 9);
        }
                      
      glEnable(GL_TEXTURE_2D);
      glActiveTexture(GL_TEXTURE6);
      glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id0[readTex]);
      shaderDoF->setLocalTexture(6);

      glBegin(GL_POLYGON);
      {
        if (normal) {
          float* v0 = &vertex[(k+0)*3];
          float* v1 = &vertex[(k+1)*3];
          float* v2 = &vertex[(k+2)*3];
          Vector dv1(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
          Vector dv2(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);
          Vector n = Cross(dv1, dv2);   n.normalize();
          glNormal3f(n.x(), n.y(), n.z());
        }

        // draw the number of vertices making up the polygon
        for(int j=0; j<poly[i]; j++) {
          float* t = &texcoord[(k+j)*3];
          float* v = &vertex[(k+j)*3];

          #if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
          #ifdef _WIN32
          if (glMultiTexCoord3f) {
          #endif // _WIN32
            glMultiTexCoord3f(GL_TEXTURE0, t[0],   t[1],     t[2]);
            if (fog) {
              float vz = mvmat[2]*v[0] + mvmat[6]*v[1] + mvmat[10]*v[2] + mvmat[14];
              glMultiTexCoord3f(GL_TEXTURE1, -vz, 0.0, 0.0);
            }
          #ifdef _WIN32
          }
          else{
            glTexCoord3f(t[0], t[1], t[2]);
          }
          #endif // _WIN32
          #else
          glTexCoord3f(t[0], t[1], t[2]);
          #endif
          glVertex3f(v[0], v[1], v[2]);
        }
      }
      glEnd();
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

      // swap read and write textures
      writeTex = 1 - writeTex;
      readTex = 1 - readTex;

      k -= poly[i];

      //writeTextureToPPM(("/home/pascal/Desktop/buffer/buffer_eye_pass1_read _" + toStr(z)).c_str(), eyeBuffer_tex_id0[readTex], textureDim_,textureDim_, GL_RGBA, GL_UNSIGNED_BYTE);
      //writeTextureToPPM(("/home/pascal/Desktop/buffer/buffer_eye_pass1_write_" + toStr(z)).c_str(), eyeBuffer_tex_id0[writeTex], textureDim_,textureDim_, GL_RGBA, GL_UNSIGNED_BYTE);
    }
    lastWritePhase1 = readTex;
    //writeTextureToPPM("/home/pascal/Desktop/lastWritePhase1" , eyeBuffer_tex_id0[readTex], textureDim_,textureDim_, GL_RGBA, GL_UNSIGNED_BYTE);

        
  // Phase 2: fp to back (front-to-back rendering)
    location = glGetUniformLocation(programId, "traversalDirection");
    glUniform1i(location, BACK_TO_FRONT);

    writeTex = 1; readTex = 0;
    z = poly.size();  

    for (unsigned int i=0, k=0; i<=zf_slicer; i++){
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frmBuffer_);
      glDrawBuffer(attachmentpoints[2+writeTex]);

      z--;
      if (dof_focusMode_ == DOF_USER){
        if ( (z >= (zf-focusRange)) && (z <= (zf+focusRange)) )
          radS = 1.0;
        else
          radS = 4.0;
      }
      else{
        Cz = dof_aperture_ * ((z-zf)/(float)z);
        if (radS <= dof_threshold_)
          radS = 1.0;
        else
          radS = Cz/dof_threshold_;
      }
    
      if (radS <= 2.0){
        location = glGetUniformLocation(programId, "numSamples");
        glUniform1i(location, 1);
      }
      else
        if (radS <= 5.0){
          location = glGetUniformLocation(programId, "numSamples");
          glUniform1i(location, 5);
        }
        else{
          location = glGetUniformLocation(programId, "numSamples");
          glUniform1i(location, 9);
        }
                
        glEnable(GL_TEXTURE_2D);
        location = glGetUniformLocation(programId, "tex6");
        glActiveTexture(GL_TEXTURE6);
        glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id1[readTex]);
        glUniform1i(location, 6);

        glBegin(GL_POLYGON);
        {
          if (normal) {
            float* v0 = &vertex[(k+0)*3];
            float* v1 = &vertex[(k+1)*3];
            float* v2 = &vertex[(k+2)*3];
            Vector dv1(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
            Vector dv2(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);
            Vector n = Cross(dv1, dv2);   n.normalize();
            glNormal3f(n.x(), n.y(), n.z());
          }

          for(int j=0; j<poly[i]; j++) {
            float* t = &texcoord[(k+j)*3];
            float* v = &vertex[(k+j)*3];

            #if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
            #ifdef _WIN32
            if (glMultiTexCoord3f) {
            #endif // _WIN32
              glMultiTexCoord3f(GL_TEXTURE0, t[0], t[1], t[2]);
              if (fog) {
                float vz = mvmat[2]*v[0] + mvmat[6]*v[1] + mvmat[10]*v[2] + mvmat[14];
                glMultiTexCoord3f(GL_TEXTURE1, -vz, 0.0, 0.0);
              }
            #ifdef _WIN32
            }
            else {
              glTexCoord3f(t[0], t[1], t[2]);
            }
            #endif // _WIN32
            #else
            glTexCoord3f(t[0], t[1], t[2]);
            #endif
            glVertex3f(v[0], v[1], v[2]);
          }
        }
        glEnd();

        // swap read and write textures
        writeTex = 1 - writeTex;
        readTex = 1 - readTex;

        k += poly[i]; 

        //writeTextureToPPM(("/home/pascal/Desktop/buffer/buffer_eye_pass2_read _" + toStr(z)).c_str(), eyeBuffer_tex_id1[readTex], textureDim_,textureDim_, GL_RGBA, GL_UNSIGNED_BYTE);
        //writeTextureToPPM(("/home/pascal/Desktop/buffer/buffer_eye_pass2_write_" + toStr(z)).c_str(), eyeBuffer_tex_id1[writeTex], textureDim_,textureDim_, GL_RGBA, GL_UNSIGNED_BYTE);
    }
    lastWritePhase2 = readTex;
    //writeTextureToPPM("/home/pascal/Desktop/lastWritePhase2" , eyeBuffer_tex_id1[readTex], textureDim_,textureDim_, GL_RGBA, GL_UNSIGNED_BYTE);

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    shaderDoF->release();
        
       
      

    // Final Phase: Render to texture
    glEnable(GL_BLEND);
    glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
    glPushMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
        glLoadIdentity();
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
          glLoadIdentity();
          programId = shaderTexDisp->activate();

            glEnable(GL_TEXTURE_2D);

            location = glGetUniformLocation(programId, "tex7");
            glActiveTexture(GL_TEXTURE7);
            glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id0[lastWritePhase1]);
            glUniform1i(location, 7);

            location = glGetUniformLocation(programId, "tex8");
            glActiveTexture(GL_TEXTURE8);
            glBindTexture(GL_TEXTURE_2D, eyeBuffer_tex_id1[lastWritePhase2]);
            glUniform1i(location, 8);
            
            location = glGetUniformLocation(programId, "option");
            glUniform1i(location, 1);

            glBegin (GL_QUADS);
              glTexCoord2f(0.0, 0.0);  glVertex3f (-1, -1, 0.0);
              glTexCoord2f(1.0, 0.0);  glVertex3f ( 1, -1, 0.0);
              glTexCoord2f(1.0, 1.0);  glVertex3f ( 1,  1, 0.0);
              glTexCoord2f(0.0, 1.0);  glVertex3f (-1,  1, 0.0);
            glEnd ();

          shaderTexDisp->release();

        glPopMatrix ();
        glMatrixMode (GL_MODELVIEW);
      glPopMatrix ();
    glPopMatrix();

  glPopMatrix();
  glPopAttrib();

  CHECK_OPENGL_ERROR();
}




void
TextureRenderer::colormap2_hardware_destroy_buffers()
{
  if (cmap2_widget_tex_id_) {
    glDeleteTextures(1, &cmap2_widget_tex_id_);
    cmap2_widget_tex_id_ = 0;
  }

  if (cmap2_tex_id_) {
    glDeleteTextures(1, &cmap2_tex_id_);
    cmap2_tex_id_ = 0;
  }

  if (cmap2_widget_framebuffer_) {
    glDeleteFramebuffersEXT(1, &cmap2_widget_framebuffer_);
    cmap2_widget_framebuffer_ = 0;
  }

  if (cmap2_framebuffer_) {
    glDeleteFramebuffersEXT(1, &cmap2_framebuffer_);
    cmap2_framebuffer_ = 0;
  }
}

  

void
TextureRenderer::colormap2_hardware_rasterize_setup()
{
  CHECK_OPENGL_ERROR();
  unsigned int s = cmap2_width_ * Pow2(cmap2_.size());
  
  if (tex_width_ != s) {
    //need to regenerate colormap2 framebuffer object.
    colormap2_hardware_destroy_buffers();
  }

  tex_width_ = s;

  if (!shader_factory_) {
    shader_factory_ = new CM2ShaderFactory();
    CHECK_OPENGL_ERROR();
    if (cmap2_shader_glsl_->create()) /* true shader->create() means it
           failed (backwards) */
    {
      cmap2_shader_glsl_->destroy();
      delete shader_factory_;
      shader_factory_ = 0;
      colormap2_hardware_destroy_buffers();
      CHECK_OPENGL_ERROR();
      return;
    }
  }
  CHECK_OPENGL_ERROR();
}
  
  

void
TextureRenderer::build_colormap2()
{
  CHECK_OPENGL_ERROR();
  if (!ShaderProgramARB::shaders_supported()) return;
  if (!cmap2_dirty_ && !alpha_dirty_) return;

  colormap2_hardware_rasterize_setup();
  colormap2_hardware_rasterize();

  cmap2_dirty_ = false;
  alpha_dirty_ = false;
  CHECK_OPENGL_ERROR();
}

void
TextureRenderer::bind_colormap1(vector<float> cmap_array,
                                unsigned int cmap_tex)
{
#if defined( GL_TEXTURE_COLOR_TABLE_SGI ) && defined(__sgi)
  glEnable(GL_TEXTURE_COLOR_TABLE_SGI);
  glColorTable(GL_TEXTURE_COLOR_TABLE_SGI,
               GL_RGBA,
               256,
               GL_RGBA,
               GL_FLOAT,
               &cmap_array[0]);
#elif defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
  if (ShaderProgramARB::shaders_supported() && glActiveTexture)
  {
    glActiveTexture(GL_TEXTURE2_ARB);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D, cmap_tex);
    // enable data texture unit 1
    glActiveTexture(GL_TEXTURE1_ARB);  
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_3D);
    glActiveTexture(GL_TEXTURE0_ARB);
  }
  else
#endif
#ifndef __sgi
#  if defined(GL_EXT_shared_texture_palette) && !defined(__APPLE__)
  {
    glEnable(GL_SHARED_TEXTURE_PALETTE_EXT);
    glColorTable(GL_SHARED_TEXTURE_PALETTE_EXT,
                 GL_RGBA,
                 256,
                 GL_RGBA,
                 GL_FLOAT,
                 &cmap_array[0]);
  }
#  else
  {
    static bool warned = false;
    if( !warned ) {
      std::cerr << "No volume colormaps available." << std::endl;
      warned = true;
    }
  }
#  endif
#endif
  CHECK_OPENGL_ERROR();
}


void
TextureRenderer::bind_colormap2()
{
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
  if (ShaderProgramARB::shaders_supported() && glActiveTexture)
  {
    // bind texture to unit 2
    glActiveTexture(GL_TEXTURE3_ARB);   // 2D transnfer fn changed from GL_TEXTURE2_ARB as this is now using the new shaders
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D, cmap2_tex_id_);
    // enable data texture unit 1
    glActiveTexture(GL_TEXTURE1_ARB);   
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_3D);
    glActiveTexture(GL_TEXTURE0_ARB);
  }
#endif
}


void
TextureRenderer::release_colormap1()
{
#if defined(GL_TEXTURE_COLOR_TABLE_SGI) && defined(__sgi)
  glDisable(GL_TEXTURE_COLOR_TABLE_SGI);
#elif defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
  if (ShaderProgramARB::shaders_supported() && glActiveTexture)
  {
    // bind texture to unit 2
    glActiveTexture(GL_TEXTURE2_ARB);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glDisable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D, 0);
    // enable data texture unit 1
    glActiveTexture(GL_TEXTURE1_ARB);  
    glDisable(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, 0);
    glActiveTexture(GL_TEXTURE0_ARB);
  }
  else
#endif
#ifndef __sgi
#  if defined(GL_EXT_shared_texture_palette) && !defined(__APPLE__)
  {
    glDisable(GL_SHARED_TEXTURE_PALETTE_EXT);
  }
#  else
  {
    // Already warned in bind.  Do nothing.
  }
#  endif
#endif
  CHECK_OPENGL_ERROR();
}


void
TextureRenderer::release_colormap2()
{
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
  if (ShaderProgramARB::shaders_supported() && glActiveTexture)
  {
    glActiveTexture(GL_TEXTURE3_ARB);  // 2D transnfer fn changed from GL_TEXTURE2_ARB as this is now using the new shaders
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0); 
    glActiveTexture(GL_TEXTURE1_ARB); 
    glDisable(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, 0);
    glActiveTexture(GL_TEXTURE0_ARB);
  }
#endif
}

} // namespace SLIVR




