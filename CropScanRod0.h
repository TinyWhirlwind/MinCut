#ifndef SEGMESH_H
#define SEGMESH_H

#include<Mesh.h>
#include<Node.h>
using namespace core;
namespace barAbutment
{
	std::pair<std::shared_ptr<core::Mesh>, std::shared_ptr<core::Mesh>> cutMesh(core::Mesh* mesh, const core::Vec3& point, const core::Vec3& normal);
}
typedef std::map<TFace, double> Face_Double_Map;
typedef std::map<TFace, std::size_t> Face_Int_Map;


class CropScanRod
{
public:
	struct SDF_property_map
	{
		double sdf_value;
		bool output_cluster_ids = false;
		bool postprocess = false;
		std::size_t label = -1;
	};
	/*double smoothing_lambda = 0.26;
	std::size_t numbers_of_rays = 1;
	double cone_angle = 2.0 / 3.0 * M_PI;
	std::size_t number_of_clusters = 2;*/
	//struct SDF_property_map
	//{
	//	double sdf_value = 0.0f;//sdf无结果用-1表示
	//	Face_Int_Map face_int_map;
	//	double cone_angle = 2.0 / 3.0 * M_PI;
	//	std::size_t numbers_of_rays = 1;
	//	bool postprocess = false;
	//	std::size_t number_of_clusters = -1;
	//	double smoothing_lambda = 0.26;
	//	bool output_cluster_ids = false;
	//};

	struct Edge
	{
		TFace* to_point;
		float weight;
	};
	

	struct ProPts
	{
		Vec3 sideP;
		Vec3 midP;
		Vec3 norP;
	};

public:
	CropScanRod(SmartMesh pDrawable, int typeNum);
	~CropScanRod(void);


public:
	bool apply();
	void setNode(SmartNode node);
	void setBoxRadius(const double boxRadius);
	void setSmoothLambda(float smoothing_lambda);
	//void logNormalizeSdfValue();
	void assignSegment(std::size_t number_of_clusters, std::map<TFace*, double> segments);


private:
	class PImpl;
	shared_ptr<PImpl> impl_;

};
#endif
