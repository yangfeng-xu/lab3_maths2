#include "Matrix4x4.hpp"
#include <cmath>
#include <stdexcept>

#define TOL 1e-6
const float EPSILON = 1e-5f;

Matrix4x4 Matrix4x4::Identity()
{
    Matrix4x4 I;
    I.At(0, 0) = 1; I.At(1, 1) = 1; I.At(2, 2) = 1; I.At(3, 3) = 1;
    return I;
}

Matrix4x4 Matrix4x4::Multiply(const Matrix4x4& B) const
{
    Matrix4x4 C{};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            double sum = 0.0;
            for (int k = 0; k < 4; ++k) {
                sum += At(i, k) * B.At(k, j);
            }
            C.At(i, j) = sum;
        }
    }
    return C;
}

Vec4 Matrix4x4::Multiply(const Vec4& v) const
{
    Vec4 res;
    res.x = At(0, 0) * v.x + At(0, 1) * v.y + At(0, 2) * v.z + At(0, 3) * v.w;
    res.y = At(1, 0) * v.x + At(1, 1) * v.y + At(1, 2) * v.z + At(1, 3) * v.w;
    res.z = At(2, 0) * v.x + At(2, 1) * v.y + At(2, 2) * v.z + At(2, 3) * v.w;
    res.w = At(3, 0) * v.x + At(3, 1) * v.y + At(3, 2) * v.z + At(3, 3) * v.w;
    return res;
}

// --------------------------------------------------------------------------
// TODO LAB 3
// --------------------------------------------------------------------------

bool Matrix4x4::IsAffine() const
{
    return std::abs(At(3, 0)) < EPSILON &&
        std::abs(At(3, 1)) < EPSILON &&
        std::abs(At(3, 2)) < EPSILON &&
        std::abs(At(3, 3) - 1.0f) < EPSILON;
}

Vec3 Matrix4x4::TransformPoint(const Vec3& p) const
{
    return Vec3(
        At(0, 0) * p.x + At(0, 1) * p.y + At(0, 2) * p.z + At(0, 3),
        At(1, 0) * p.x + At(1, 1) * p.y + At(1, 2) * p.z + At(1, 3),
        At(2, 0) * p.x + At(2, 1) * p.y + At(2, 2) * p.z + At(2, 3)
    );
}

Vec3 Matrix4x4::TransformVector(const Vec3& v) const
{
    return Vec3(
        At(0, 0) * v.x + At(0, 1) * v.y + At(0, 2) * v.z,
        At(1, 0) * v.x + At(1, 1) * v.y + At(1, 2) * v.z,
        At(2, 0) * v.x + At(2, 1) * v.y + At(2, 2) * v.z
    );
}

Matrix4x4 Matrix4x4::Translate(const Vec3& t)
{
 
    Matrix4x4 M = Identity();
    M.At(0, 3) = t.x;
    M.At(1, 3) = t.y;
    M.At(2, 3) = t.z;
    return M;
}

Matrix4x4 Matrix4x4::Scale(const Vec3& s)
{
    Matrix4x4 M = Identity();
    M.At(0, 0) = s.x;
    M.At(1, 1) = s.y;
    M.At(2, 2) = s.z;
    return M;
}

Matrix4x4 Matrix4x4::Rotate(const Matrix3x3& R)
{
    Matrix4x4 M = Identity();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            // Assuming Matrix3x3 uses At() or similar access
            M.At(i, j) = R.At(i, j);
        }
    }
    return M;
}

Matrix4x4 Matrix4x4::Rotate(const Quat& q)
{
    Matrix4x4 M = Identity();

    // Pre-calculate repeated values to save performance
    float xx = q.x * q.x;
    float yy = q.y * q.y;
    float zz = q.z * q.z;
    float xy = q.x * q.y;
    float xz = q.x * q.z;
    float yz = q.y * q.z;
    float wx = q.s * q.x;
    float wy = q.s * q.y;
    float wz = q.s * q.z;

    // Fill the top-left 3x3 rotation matrix
    M.At(0, 0) = 1.0f - 2.0f * (yy + zz);
    M.At(0, 1) = 2.0f * (xy - wz);
    M.At(0, 2) = 2.0f * (xz + wy);

    M.At(1, 0) = 2.0f * (xy + wz);
    M.At(1, 1) = 1.0f - 2.0f * (xx + zz);
    M.At(1, 2) = 2.0f * (yz - wx);

    M.At(2, 0) = 2.0f * (xz - wy);
    M.At(2, 1) = 2.0f * (yz + wx);
    M.At(2, 2) = 1.0f - 2.0f * (xx + yy);

    return M;
}

Matrix4x4 Matrix4x4::FromTRS(const Vec3& t, const Matrix3x3& R, const Vec3& s)
{
    Matrix4x4 M = Identity();

    // 1. Apply Rotation and Scale to the top-left 3x3
    for (int r = 0; r < 3; ++r) {
        M.At(r, 0) = R.At(r, 0) * s.x;
        M.At(r, 1) = R.At(r, 1) * s.y;
        M.At(r, 2) = R.At(r, 2) * s.z;
    }

    // 2. Apply Translation to the last column
    M.At(0, 3) = t.x;
    M.At(1, 3) = t.y;
    M.At(2, 3) = t.z;

    return M;
}

Matrix4x4 Matrix4x4::FromTRS(const Vec3& t, const Quat& q, const Vec3& s)
{
    Matrix4x4 M = Identity();

    // 1. Calculate Quaternion terms
    float xx = q.x * q.x;
    float yy = q.y * q.y;
    float zz = q.z * q.z;
    float xy = q.x * q.y;
    float xz = q.x * q.z;
    float yz = q.y * q.z;
    float wx = q.s * q.x;
    float wy = q.s * q.y;
    float wz = q.s * q.z;

    // 2. Apply Rotation * Scale
    // Column 0
    M.At(0, 0) = (1.0f - 2.0f * (yy + zz)) * s.x;
    M.At(1, 0) = (2.0f * (xy + wz)) * s.x;
    M.At(2, 0) = (2.0f * (xz - wy)) * s.x;

    // Column 1
    M.At(0, 1) = (2.0f * (xy - wz)) * s.y;
    M.At(1, 1) = (1.0f - 2.0f * (xx + zz)) * s.y;
    M.At(2, 1) = (2.0f * (yz + wx)) * s.y;

    // Column 2
    M.At(0, 2) = (2.0f * (xz + wy)) * s.z;
    M.At(1, 2) = (2.0f * (yz - wx)) * s.z;
    M.At(2, 2) = (1.0f - 2.0f * (xx + yy)) * s.z;

    // 3. Apply Translation (Column 3)
    M.At(0, 3) = t.x;
    M.At(1, 3) = t.y;
    M.At(2, 3) = t.z;

    return M;
}

Matrix4x4 Matrix4x4::InverseTR() const
{
    Matrix4x4 M;
    return M;
}

Matrix4x4 Matrix4x4::InverseTRS() const
{
    Matrix4x4 M;
    return M;
}

Vec3 Matrix4x4::GetTranslation() const
{
    return Vec3{};
}

Matrix3x3 Matrix4x4::GetRotationScale() const
{
    Matrix3x3 M;
    return M;
}

Vec3 Matrix4x4::GetScale() const
{
	return Vec3{};
}

Matrix3x3 Matrix4x4::GetRotation() const
{
    Matrix3x3 M;
    return M;
}

Quat Matrix4x4::GetRotationQuat() const
{
	return Quat{};
}

void Matrix4x4::SetTranslation(const Vec3& t)
{
    
}

void Matrix4x4::SetScale(const Vec3& s)
{
    
}

void Matrix4x4::SetRotation(const Matrix3x3& R)
{
    
}

void Matrix4x4::SetRotation(const Quat& q)
{
   
}

void Matrix4x4::SetRotationScale(const Matrix3x3& RS)
{
    
}