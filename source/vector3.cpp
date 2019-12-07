#include "vector3.h"

Vec3d::~Vec3d(void)
{

}

Vec3d operator - ( const Vec3d& A){
	return Vec3d( -A.x, -A.y , -A.z );
}

Vec3d operator + ( const Vec3d& A , const Vec3d& B ){
	return Vec3d( A.x + B.x , A.y + B.y , A.z + B.z );
}

Vec3d operator - ( const Vec3d& A , const Vec3d& B ){
	return Vec3d( A.x - B.x , A.y - B.y , A.z - B.z );
}

Vec3d operator * ( const Vec3d& A , const Vec3d& B ){
	return Vec3d( A.x * B.x , A.y * B.y , A.z * B.z );
}

Vec3d operator / ( const Vec3d& A , const Vec3d& B ){
	return Vec3d( A.x / B.x , A.y / B.y , A.z / B.z );
}

Vec3d operator + ( const Vec3d& A , const double& t ){
	return Vec3d( A.x + t , A.y + t , A.z + t );
}

Vec3d operator - ( const Vec3d& A , const double& t ){
	return Vec3d( A.x - t , A.y - t , A.z - t );
}

Vec3d operator * ( const Vec3d& A , const double& t ){
	return Vec3d( A.x * t , A.y * t , A.z * t );
}

Vec3d operator / ( const Vec3d& A , const double& t ){
	return Vec3d( A.x / t , A.y / t , A.z / t );
}

Vec3d Vec3d::getUnitVectorOfThis(){
	double module = sqrt(x * x + y * y + z * z);
	if(module < Config::EPS)
		return Vec3d(0,0,1);
	return *this / module;
}

Vec3d Vec3d::getVerticalVectorCrossed001(){
	//tmp is vertical to this and (0,0,1)
	Vec3d tmp  = Vec3d(y,- x,0 );
	double eps = Config::EPS;
	if(fabs(x) < eps && fabs(y) < eps && fabs(z) < eps)
		return Vec3d(1,0,0);
	return tmp.getUnitVectorOfThis();
	
}

Vec3d Vec3d::getCross(const Vec3d& A){
	return Vec3d(y * A.z - z * A.y , z * A.x - x * A.z , x * A.y - y * A.x );
}

double Vec3d::getDot(const Vec3d& A){
	return x * A.x + y * A.y + z * A.z;
}





