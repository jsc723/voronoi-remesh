#ifndef MYVEC3D_H
#define MYVEC3D_H
#include "config.h"
#include <cmath>
class Vec3d
{
	
public:
	Vec3d(double X = 0 , double Y = 0 , double Z = 0 ) : x( X ) , y( Y ) , z( Z ) {};
	~Vec3d(void);
	double x,y,z;
	friend Vec3d operator - ( const Vec3d& );//ȡ��
	friend Vec3d operator + ( const Vec3d& , const Vec3d& );
	friend Vec3d operator - ( const Vec3d& , const Vec3d& );
	friend Vec3d operator * ( const Vec3d& , const Vec3d& );
	friend Vec3d operator / ( const Vec3d& , const Vec3d& );
	friend Vec3d operator + ( const Vec3d& , const double& );
	friend Vec3d operator - ( const Vec3d& , const double& );
	friend Vec3d operator * ( const Vec3d& , const double& );
	friend Vec3d operator / ( const Vec3d& , const double& );

	Vec3d getUnitVectorOfThis();
	Vec3d getVerticalVectorCrossed001();
	Vec3d getCross(const Vec3d&);
	double getDot(const Vec3d&);
};

#endif

