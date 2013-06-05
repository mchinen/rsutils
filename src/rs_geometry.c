/*
 *  rs_geometry.c
 *  RS
 *
 *  Created by Administrator on 12/21/10.
 *  Copyright 2010 Roughsoft. All rights reserved.
 *
 */

#include <math.h>
#include <float.h>

#include "rs_util.h"

#include "rs_geometry.h"

void rs_overlap_result_flip(int *inoutResult)
{
    int newResult = 0;
    
    if (*inoutResult & eRSCubeOverlapResultTop)
        newResult |= eRSCubeOverlapResultBottom;
    if (*inoutResult & eRSCubeOverlapResultBottom)
        newResult |= eRSCubeOverlapResultTop;
    if (*inoutResult & eRSCubeOverlapResultRight)
        newResult |= eRSCubeOverlapResultLeft;
    if (*inoutResult & eRSCubeOverlapResultLeft)
        newResult |= eRSCubeOverlapResultRight;
    *inoutResult = newResult;
}

void rs_overlap_result_normal(int result, RSVec3* outNormal)
{
   rs_vec_zero(outNormal);
   switch (result) {
   case eRSCubeOverlapResultTop:
      outNormal->y = 1.;
      break;
   case eRSCubeOverlapResultBottom:
      outNormal->y = -1.;
      break;
   case eRSCubeOverlapResultRight:
      outNormal->x = 1.;
      break;
   case eRSCubeOverlapResultLeft:
      outNormal->x = -1.;
      break;
   default:
      break;
   }
}

int rs_cube_overlaps(RSCube *a, RSCube *b, int *outResult)
{
    int res;
    RSVec3 atob;
    rs_vec_subtract(&b->pos, &a->pos, &atob);
    
    *outResult = eRSCubeOverlapResultNone;
    res =   fabs(atob.x) <= (a->size.x / 2 + b->size.x / 2) &&
            fabs(atob.y) <= (a->size.y / 2 + b->size.y / 2) &&
            fabs(atob.z) <= (a->size.z / 2 + b->size.z / 2);
            
    if (res) {
        //hit the right side of A
        *outResult = atob.x > 0 ? eRSCubeOverlapResultRight : eRSCubeOverlapResultLeft;
        // now note how far in each axis is
        // we expect the actual collision side to be the *least* overlapping
        if ((a->size.x / 2 + b->size.x / 2) - fabs(atob.x) >
            (a->size.y / 2 + b->size.y / 2) - fabs(atob.y)) {
            *outResult = atob.y < 0 ? eRSCubeOverlapResultBottom : eRSCubeOverlapResultTop;
        }

    }
    return res;
}

void rs_vec_zero(RSVec3 *a)
{
    a->x = a->y = a->z = 0.0f;
}

void rs_vec_scale(RSVec3 *a, float scalar, RSVec3 *out)
{
    out->x = a->x * scalar;
    out->y = a->y * scalar;
    out->z = a->z * scalar;    
}

void rs_vec_add(RSVec3 *a, RSVec3 *b, RSVec3 *out)
{
    out->x = a->x + b->x;
    out->y = a->y + b->y;
    out->z = a->z + b->z;    
}

void rs_vec_subtract(RSVec3 *a, RSVec3 *b, RSVec3 *out)
{
    out->x = a->x - b->x;
    out->y = a->y - b->y;
    out->z = a->z - b->z;
}
//euclidean distance.
float rs_vec_dist(RSVec3 *a, RSVec3 *b)
{
    RSVec3 d;
    rs_vec_subtract(a, b, &d);
    return sqrt((d.x * d.x) + (d.y * d.y) + (d.z * d.z));
}

float rs_vec_length(RSVec3 *a)
{
    RSVec3 zero = {0.0f, 0.0f, 0.0f};
    return rs_vec_dist(a, &zero);
}

float rs_vec_dot(RSVec3 *a, RSVec3 *b)
{
   return a->x * b->x + a->y * b->y + a->z * b->z;
}

void rs_vec_cross(RSVec3 *a, RSVec3 *b, RSVec3 *out)
{
   out->x = a->y * b->z - b->y * a->z;
   out->y = a->z * b->x - b->z * a->x;
   out->z = a->x * b->y - b->x * a->y;
}

float rs_vec_norm(RSVec3 *in, RSVec3 *out)
{
    float dist = rs_vec_length(in);
    *out = *in;
    if (dist) {
        out->x /= dist;
        out->y /= dist;
        out->z /= dist;
    }
    return dist;
}

void rs_vec_abs(RSVec3 *in, RSVec3 *out)
{
   out->x = fabs(in->x);
   out->y = fabs(in->y);
   out->z = fabs(in->z);
}

void rs_vec_neg(RSVec3 *in, RSVec3 *out)
{
   out->x = -1.0 * in->x;
   out->y = -1.0 * in->y;
   out->z = -1.0 * in->z;
}

//rotates on xy plane (z stays same)
void rs_vec_rotatexy(RSVec3 *in, float radAngle, RSVec3 *out)
{
   // cache in case in == out
   float inx = in->x;
   float cosf = cos(radAngle);
   float sinf = sin(radAngle);
   out->x = inx * cosf - in->y * sinf;
   out->y = inx * sinf + in->y * cosf;
   out->z = in->z;
}

// angle between two vectors
float rs_vec_anglexy(RSVec3 *a, RSVec3 *b, RSVec3 *ref)
{
   /*
   float dx21 = a->x - ref->x;
   float dx31 = b->x - ref->x;
   float dy21 = a->y - ref->y;
   float dy31 = b->y - ref->y;
   float m12 = sqrt( dx21*dx21 + dy21*dy21 );
   float m13 = sqrt( dx31*dx31 + dy31*dy31 );

   if (m12 * m13 == 0)
      return 0.;

   return acos( (dx21*dx31 + dy21*dy31) / (m12 * m13) );
*/

   RSVec3 aNorm, bNorm;

         
   rs_vec_subtract(a, ref, &aNorm);
   rs_vec_subtract(b, ref, &bNorm);

   float val = atan2(bNorm.y, bNorm.x) - atan2(aNorm.y, aNorm.x);
   if (val > M_PI)
      val -= 2 * M_PI;
   else if (val < -M_PI)
      val += 2* M_PI;
   return val;

   /*

   // assume no div by zero input
   rs_vec_norm(&aNorm, &aNorm);
   rs_vec_norm(&bNorm, &bNorm);
   rs_vec_cross(a, b, &cross);
   return (atan2(rs_vec_length(&cross), rs_vec_dot(&aNorm, &bNorm)) *
           (cross.z > 0 ? 1 : -1));
   */
}

/*
void rs_vec_rotatexy_offcenter(RSVec3 *in, float radAngle, RSVec3 *out)
{
   
}
*/

 // only works for XYZ axis aligned cube
void rs_cube_rotatexy_bounds(RSCube* in, float radangle, RSCube* out)
{
   // reconstruct the 4 rect points
   RSVec3 points[4];
   points[0].x = points[3].x = -in->size.x / 2;
   points[1].x = points[2].x =  in->size.x / 2;
   
   points[0].y = points[1].y = -in->size.y / 2;
   points[2].y = points[3].y =  in->size.y / 2;

   points[1].z = points[2].z = points[3].z = points[0].z = in->pos.z;

   // rotate each point
   float xmin, xmax, ymin, ymax;
   xmin = ymin = FLT_MAX;
   xmax = ymax = -FLT_MAX;

   unsigned int i;
   for (i = 0; i < 4; i++) {
      rs_vec_rotatexy(&points[i], radangle, &points[i]);
      xmin = SIMP_MIN(points[i].x, xmin);
      xmax = SIMP_MAX(points[i].x, xmax);
      ymin = SIMP_MIN(points[i].y, ymin);
      ymax = SIMP_MAX(points[i].y, ymax);
   }

   // construct a new cube using the min/max of x and y
   out->pos = in->pos;
   out->size.x = xmax - xmin;
   out->size.y = ymax - ymin;
   out->size.z = in->size.z;
}

void rs_cube_rotated_points(RSCube* in, float radangle, RSVec3 points[4])
{
   points[0].x = points[3].x = -in->size.x / 2;
   points[1].x = points[2].x = in->size.x / 2;
   
   points[0].y = points[1].y = -in->size.y / 2;
   points[2].y = points[3].y = in->size.y / 2;   

   points[1].z = points[2].z = points[3].z = points[0].z = in->pos.z;
   unsigned int i;
   for (i = 0; i < 4; i++) {
      rs_vec_rotatexy(&points[i], radangle, &points[i]);
      points[i].x += in->pos.x;
      points[i].y += in->pos.y;
   }
}

float rs_overlap_result_radians(int result)
{
   switch (result) {
   case eRSCubeOverlapResultRight:
      return M_PI * 3 / 2;
   case eRSCubeOverlapResultBottom:
      return M_PI;
   case eRSCubeOverlapResultLeft:
      return M_PI / 2;
   case eRSCubeOverlapResultTop:
   default:
      return 0.;
   }
}

float rs_vec_project_axis_vec(RSVec3 *in, RSVec3* axis)
{
   float temp = ((in->x * axis->x + in->y * axis->y)
                 / (axis->x * axis->x + axis->y * axis->y));
   // add dot product of both
   return temp * axis->x * axis->x + temp * axis->y * axis->y;
}

int rs_cube_overlaps_rotated(RSCube *a, float aAngleRad,
                                 RSCube *b, float bAngleRad,
                                 int *outResult)
{
   RSVec3 pointsA[4], pointsB[4];
   RSVec3 axis[4];
   float min[2], max[2];
   float proj, temp;
   float centerA, centerB, sizeA, sizeB;
   RSVec3* points;
   float minOverlap = FLT_MAX;
   unsigned int i,j, k;

   rs_cube_rotated_points(a, aAngleRad, pointsA);
   rs_cube_rotated_points(b, bAngleRad, pointsB);
   
   // get the two axis for a.
   // we take the bottom axis and then the left axis
   rs_vec_subtract(&pointsA[1], &pointsA[0], &axis[0]);
   rs_vec_subtract(&pointsA[3], &pointsA[0], &axis[1]);

   rs_vec_subtract(&pointsB[1], &pointsB[0], &axis[2]);
   rs_vec_subtract(&pointsB[3], &pointsB[0], &axis[3]);

   // normalizing the axis vectors is not needed for collision detection
   // but it can make the distance comparisons done at the end of the function
   // more meaningful to determine which side we have hit
   rs_vec_norm(&axis[0], &axis[0]);
   rs_vec_norm(&axis[1], &axis[1]);
   rs_vec_norm(&axis[2], &axis[2]);
   rs_vec_norm(&axis[3], &axis[3]);

   // go over each of the axis
   // 1st min is for A, 2nd is for B
   for (i = 0; i < 4; i++) {
      min[0] = min[1] = FLT_MAX;
      max[0] = max[1] = -FLT_MAX;
      // for each axis, go over each point
      for (j = 0; j < 4; j++) {
         for (k = 0; k < 2; k++) {
            points = k ? pointsB : pointsA;

            // do the projection of one corner point to the current axis
            proj = rs_vec_project_axis_vec(&points[j], &axis[i]);
            min[k] = SIMP_MIN(min[k], proj);
            max[k] = SIMP_MAX(max[k], proj);
         }
      }
      // now see if we have a 1 dimensional overlap.
      // if we don't there is no overlap and we should bail
      // use the same center-based overlap detection as rs_cube_overlaps
      centerA = (min[0] + max[0]) / 2;
      centerB = (min[1] + max[1]) / 2;
      sizeA   = max[0] - min[0];
      sizeB   = max[1] - min[1];
      temp = centerB - centerA;
      if (!(fabs(temp) <= sizeA / 2 + sizeB / 2)) {
         // if we didn't overlap there is no collsion so bail
         *outResult = eRSCubeOverlapResultNone;
         return 0;
      } else if (i < 2) {
         // see if we have a collision.  orientation is from the perspective
         // of A, so we only do this for the 1st and 2nd axis checks
         // figure out how much we are overlapping
         float overlapMax, overlapMin, overlap;
         overlapMax = SIMP_MIN(max[0], max[1]);
         overlapMin = SIMP_MAX(min[0], min[1]);
         overlap = overlapMax - overlapMin;

         if (fabs(overlap) < minOverlap) {
            if (i == 0) {
               *outResult = temp > 0 ? eRSCubeOverlapResultRight : eRSCubeOverlapResultLeft;
            } else {
               *outResult = temp < 0 ? eRSCubeOverlapResultBottom : eRSCubeOverlapResultTop;
            }
            minOverlap = fabs(overlap);
         }
      }
   }

   return 1;//TODO:replace
}
