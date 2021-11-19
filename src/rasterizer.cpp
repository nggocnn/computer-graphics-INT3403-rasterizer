#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    rgb_framebuffer_target[3 * (y * width + x)    ] = (unsigned char)(c.r * 255);
    rgb_framebuffer_target[3 * (y * width + x) + 1] = (unsigned char)(c.g * 255);
    rgb_framebuffer_target[3 * (y * width + x) + 2] = (unsigned char)(c.b * 255);
  }

  void RasterizerImp::fill_supersample(size_t x, size_t y, size_t s, Color c) {
    // TODO: Task 2: You may want to implement this function. Hint: our solution uses one line
    auto sx = (int)floor(x);
    auto sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;
    sample_buffer[sample_rate * (sy * width + sx) + s] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    auto sx = (int)floor(x);
    auto sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    for (int s = 0; s < sample_rate; s++) {
      fill_supersample(sx, sy, s, color);
    }
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
                                     float x1, float y1,
                                     Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
                                         float x1, float y1,
                                         float x2, float y2,
                                         Color color) {
    // Create z_hat vector for convenience in computing the normal vector
    Vector3D z(0, 0, 1);
    // Create three triangle vertices
    Vector3D p0(x0, y0, 0);
    Vector3D p1(x1, y1, 0);
    Vector3D p2(x2, y2, 0);
    // Ensure that the vertices are enumerated in a clockwise fashion
    //       p0
    //      / \
    //    p2---p1
    if (cross(((p1 + p2) / 2) - p0, p1 - p0).z < 0) {
      swap(p1, p2);
    }
    // Create all three triangle edge lines
    Vector3D lin0 = p0 - p1;
    Vector3D lin1 = p1 - p2;
    Vector3D lin2 = p2 - p0;
    // Find clockwise normals for each line (these point INTO the triangle)
    Vector3D n0 = cross(z,lin0);
    Vector3D n1 = cross(z,lin1);
    Vector3D n2 = cross(z,lin2);
    // Bounding box limits for triangle
    float x_list[] = { x0,x1,x2 };
    float y_list[] = { y0,y1,y2 };
    // x pixel bounds
    int min_x = (int)floor(*std::min_element(x_list, x_list+3));
    int max_x = (int)ceil(*std::max_element(x_list, x_list+3));
    // y pixel bounds
    int min_y = (int)floor(*std::min_element(y_list, y_list+3));
    int max_y = (int)ceil(*std::max_element(y_list, y_list+3));

    // Get sample rate
    int rate = (int)sqrt(sample_rate);   // sample rate in 1D

    // Loop through pixels in bounding box of triangle
    for (int y = min_y; y < max_y; y++) {
      for (int x = min_x; x < max_x; x++) {
        Vector3D p(x, y, 0);    // screen point location
        int s = 0;              // linear index for supersample box
        for (int j = 0; j < rate; j++) {
          p.y = (float)y + ((float)j + 0.5) / (float)rate;
          for (int i = 0; i < rate; i++) {
            p.x = (float)x + ((float)i + 0.5) / (float)rate;
            // Line test: Is pixel inside or outside of the triangle?
            if ((dot(p - p1, n0) >= 0) && (dot(p - p2, n1) >= 0) && (dot(p - p0, n2) >= 0)) {
              // fill the super sampling buffer
              fill_supersample(x, y, s, color);
            }
            s++;
          }
        }
      }
    }
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling

    // TODO: Task 2: Update to implement super-sampled rasterization

  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
                                                            float x1, float y1, Color c1,
                                                            float x2, float y2, Color c2) {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
    // Create z_hat vector for convenience in computing the normal vector
    Vector3D z(0, 0, 1);
    // Create three triangle vertices
    Vector3D p0(x0, y0, 0);
    Vector3D p1(x1, y1, 0);
    Vector3D p2(x2, y2, 0);
    // Ensure that the vertices are enumerated in a clockwise fashion
    //       p0
    //      / \
    //    p2---p1
    if (cross(((p1 + p2) / 2) - p0, p1 - p0).z < 0) {
      swap(p1, p2);
    }
    // Create all three triangle edge lines
    Vector3D lin0 = p0 - p1;
    Vector3D lin1 = p1 - p2;
    Vector3D lin2 = p2 - p0;
    // Find clockwise normals for each line (these point INTO the triangle)
    Vector3D n0 = cross(z, lin0);
    Vector3D n1 = cross(z, lin1);
    Vector3D n2 = cross(z, lin2);
    // Bounding box limits for triangle
    float x_list[] = { x0,x1,x2 };
    float y_list[] = { y0,y1,y2 };
    // x pixel bounds
    int min_x = (int)floor(*std::min_element(x_list, x_list + 3));
    int max_x = (int)ceil(*std::max_element(x_list, x_list + 3));
    // y pixel bounds
    int min_y = (int)floor(*std::min_element(y_list, y_list + 3));
    int max_y = (int)ceil(*std::max_element(y_list, y_list + 3));

    // Get sample rate
    int sr = (int)sqrt(sample_rate);   // sample rate in 1D

    // Matrix for computing berrycentric coefficients
    Matrix3x3 M(x0, x1, x2, y0, y1, y2, 1, 1, 1);
    M = M.inv();

    // Loop through pixels in bounding box of triangle
    for (int y = min_y; y < max_y; y++) {
      for (int x = min_x; x < max_x; x++) {
        Vector3D p(x, y, 1);    // screen point location
        int s = 0;              // linear index for supersample box
        for (int j = 0; j < sr; j++) {
          p.y = (float)y + ((float)j + 0.5) / (float)sr;
          for (int i = 0; i < sr; i++) {
            p.x = (float)x + ((float)i + 0.5) / (float)sr;
            // Line test: Is pixel inside or outside of the triangle?
            if ((dot(p - p1, n0) >= 0) && (dot(p - p2, n1) >= 0) && (dot(p - p0, n2) >= 0)) {
              // Get berrycentric coefficients
              Vector3D weights = M*p;
              weights.z = 1 - weights.x - weights.y;
              // fill the super sampling buffer
              fill_supersample(x, y, s, weights.x*c0 + weights.y*c1 + weights.z*c2);
            }
            s++;
          }
        }
      }
    }
  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex) {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
    // Create z_hat vector for convenience in computing the normal vector
    Vector3D z(0, 0, 1);
    // Create three screenspace triangle vertices
    Vector3D p0(x0, y0, 0);
    Vector3D p1(x1, y1, 0);
    Vector3D p2(x2, y2, 0);
    // Ensure that the vertices are enumerated in a clockwise fashion
    //       p0
    //      / \
    //    p2---p1
    if (cross(((p1 + p2) / 2) - p0, p1 - p0).z < 0) {
      swap(p1, p2);
    }

    // Create all three triangle edge lines
    Vector3D lin0 = p0 - p1;
    Vector3D lin1 = p1 - p2;
    Vector3D lin2 = p2 - p0;
    // Find clockwise normals for each line (these point INTO the triangle)
    Vector3D n0 = cross(z, lin0);
    Vector3D n1 = cross(z, lin1);
    Vector3D n2 = cross(z, lin2);
    // Bounding box limits for triangle
    float x_list[] = { x0,x1,x2 };
    float y_list[] = { y0,y1,y2 };
    // x pixel bounds
    int min_x = (int)floor(*std::min_element(x_list, x_list + 3));
    int max_x = (int)ceil(*std::max_element(x_list, x_list + 3));
    // y pixel bounds
    int min_y = (int)floor(*std::min_element(y_list, y_list + 3));
    int max_y = (int)ceil(*std::max_element(y_list, y_list + 3));

    // Get sample rate
    int rate = sqrt(sample_rate);   // sample rate in 1D

    // Instantiate SampleParams structure for setting up sampling method
    SampleParams S;
    S.lsm = lsm;
    S.psm = psm;
    // Matrix and vectors for computing berrycentric weights
    Matrix3x3 M(x0, x1, x2, y0, y1, y2, 1, 1, 1);
    M = M.inv();
    Vector3D u = Vector3D(u0, u1, u2);
    Vector3D v = Vector3D(v0, v1, v2);

    // Loop through pixels in bounding box of triangle
    for (int y = min_y; y < max_y; y++) {
      for (int x = min_x; x < max_x; x++) {
        Vector3D p(x, y, 1);    // screen point location
        int s = 0;              // linear index for supersample box
        for (int j = 0; j < rate; j++) {
          p.y = (float)y + ((float)j + 0.5) / (float)rate;
          for (int i = 0; i < rate; i++) {
            p.x = (float)x + ((float)i + 0.5) / (float)rate;

            // Line test: Is pixel inside or outside the triangle?
            if ((dot(p - p1, n0) >= 0) && (dot(p - p2, n1) >= 0) && (dot(p - p0, n2) >= 0)) {
              Vector3D w_p = M * p;                             // Get Alpha, Beta, Gamma for barycentric coords of p = (x,y)
              Vector3D w_p_dx = M * (p + Vector3D(1, 0, 0));    // Get Alpha, Beta, Gamma for barycentric coords of p_dx = (x+1,y)
              Vector3D w_p_dy = M * (p + Vector3D(0, 1, 0));    // Get Alpha, Beta, Gamma for barycentric coords of p_dy = (x, y+1)
              // Ensure alpha, beta, gamma sum to 1
              w_p.z = 1 - w_p.x - w_p.y;
              w_p_dx.z = 1 - w_p_dx.x - w_p_dx.y;
              w_p_dy.z = 1 - w_p_dy.x - w_p_dy.y;
              // Get UV space coordinates using barycentric coords
              // Set the target uv coordinate in the sample parameter object
              Vector2D p_uv = Vector2D(dot(w_p, u), dot(w_p, v));
              Vector2D p_uv_dx = Vector2D(dot(w_p_dx, u), dot(w_p_dx, v));
              Vector2D p_uv_dy = Vector2D(dot(w_p_dy, u), dot(w_p_dy, v));
              // Make differential by subtracting target uv coordinate from uv cordinates resulting mapping to (x+1,y) and (x,y+1)
              S.p_uv = p_uv;
              S.p_dx_uv = p_uv_dx - p_uv;
              S.p_dy_uv = p_uv_dy - p_uv;
              Color c = tex.sample(S);
//              Color c;
//              if (psm == P_NEAREST) {
//                c = tex.sample_nearest(p_uv);
//              } else if (psm == P_LINEAR) {
//                c = tex.sample_bilinear(p_uv);
//              }
              fill_supersample(x, y, s, c);
            }
            s++;
          }
        }
      }
    }



  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support
    this->sample_rate = rate;
    sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height) {
    // TODO: Task 2: You may want to update this function for supersampling support
    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;
    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        Color avg_color = Color::Black;
        for (int s = 0; s < sample_rate; s++) {
          avg_color += sample_buffer[sample_rate * (y * width + x) + s];
        }
        avg_color.r /= sample_rate;
        avg_color.g /= sample_rate;
        avg_color.b /= sample_rate;
        fill_pixel(x, y, avg_color);
      }

    }
  }

  Rasterizer::~Rasterizer() { }


}// CGL
