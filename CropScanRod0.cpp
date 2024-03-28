#include<CropScanRod.h>
#include <ray3.h>
#include<SmartAnn.h>
#include <MeshQuery.h>
#include <ObbFilter.h>
#include <PlaneSection.h>
#include<PointGeometry.h>
#include<LineGeometry.h>
#include<minCut.h>
#include<SelectFilter.h>
#include <chrono>
#include<OcclusionModel.h>
#include<BasisGeometryOperation.h>
#include<CurvatureFilter.h>
#include<CleanFilter.h>
/*
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/mesh_segmentation.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SM;
typedef boost::graph_traits<SM>::face_descriptor face_descriptor;*/
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/Polyhedron_items_with_id_3.h>
//#include <CGAL/mesh_segmentation.h>
//#include <CGAL/property_map.h>
//#include <iostream>
//#include <fstream>
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;
//typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
//typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
//typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/mesh_segmentation.h>
#include <SurfaceMeshProcessing.h>
#include <MeshCurvature2.h>
#include <FillHole.h>
#include <iostream>
#include <fstream>
#include <sstream>
//typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SM;
typedef boost::graph_traits<SM>::face_descriptor face_descriptor;

using namespace core;
class CropScanRod::PImpl
{
public:
	PImpl(SmartMesh pDrawable, int typeNum)
	{
		pMesh_ = pDrawable;
		typeNum_ = typeNum;
		boxRadius_ = 0.1;
		meshList_.clear();
		sdf_map.clear();
		cluster_num_ = 2; 
		numbers_of_rays_ = 25;
		ann_ = pMesh_->createTree(ITree::TT_ANN)->asAnn();
		mq_ = pMesh_->createTree(ITree::TT_QUERY)->asQuery();
		queryRadius_ = 5.f;
		scanCrodHeight_ = 14.f;
		scanRodNum_ = 6;
		/*queryRadius_ = 4.2f;
		scanCrodHeight_ = 8.f;
		scanRodNum_ = 8;*/
	}
	~PImpl() {
	};

public:
	SmartMesh pMesh_;
	SmartNode node_;
	int typeNum_;
	double boxRadius_;
	float threshold_;//阈值
	float smoothing_lambda_;//平滑系数
	std::vector<SmartMesh> meshList_;

	MeshQuery* mq_;
	MeshQuery::Intersection out_intersection_;
	Ray3d ray_;
	Vec3 planeNor_;
	Vec3 planeCenter_;
	Vec3 planeMidNor_;
	Vec3 planeSideNor_;
	SmartAnn* ann_;
	shared_ptr<ObbFilter> obbObject_;
	//================================================================
	float queryRadius_;//4.5
	float scanCrodHeight_;
	int scanRodNum_;//6
	//================================================================




	std::map<TFace*, SDF_property_map> sdf_map;
	int cluster_num_;//分类数
	int numbers_of_rays_;//射线数量
public:
	bool apply();
	void calcOcclusalSurface();
	static bool orderValue(std::pair<TVertex*, float> a, std::pair<TVertex*, float> b);
	bool bIntersection(Vec3 rayOrigin, Vec3 rayNormal);
	std::set<TVertex*> findSeedPts();
	void tracerseCrop(std::set<TVertex*> seedPts);
	void tracerseCrop2(std::set<TVertex*>& seedPts);

	void sdf_values(const double boxRadius);
	void BilateralSdfFilter(SmartMesh mesh);//双边滤波
	void XOYOZ(Vec3& zDir, Vec3& xDir, Vec3& yDir);
	void SoftCluster();//软聚类
	float DistributionProbability(std::size_t cluster);//分布概率,EM
	void EnergyFunction(float threshold, double smoothing_lambda);
	double alpha_expansion_graphcut(float threshold, double smoothing_lambda);
};
CropScanRod::CropScanRod(SmartMesh pDrawable, int typeNum)
{
	impl_.reset(new PImpl(pDrawable, typeNum));
}

CropScanRod::~CropScanRod()
{
}

bool CropScanRod::apply()
{
	return impl_->apply();
}

void CropScanRod::setNode(SmartNode node)
{
	impl_->node_ = node;
}

void CropScanRod::setBoxRadius(const double boxRadius)
{
	impl_->boxRadius_ = boxRadius;
}

void CropScanRod::PImpl::calcOcclusalSurface()
{
	Vec3f occlusalDir;
	Vec3f buccalDir;
	Vec3f archMeshCenter;
	Vec3f planeCenter;
	OcclusionModel occSurface(nullptr);
	occSurface.calcOcclusionModelDirections(pMesh_.get(), occlusalDir, buccalDir, archMeshCenter);
	occSurface.calcOcclusalPlane(pMesh_.get(), occlusalDir, buccalDir, planeCenter, OcclusalPlaneType::LOOSE);
	planeCenter_ = planeCenter;
	planeNor_ = occlusalDir;
	planeMidNor_ = buccalDir;

	ray_.SetOrigin(Vec3d::Construct(planeCenter_));
	ray_.SetDirection(Vec3d::Construct(planeNor_));
	if (mq_->intersects(ray_, out_intersection_))
	{
		planeNor_ = -planeNor_;
	}
	else
	{
		Vec3 meshCenterDir = pMesh_->asMesh()->getMeshCenter() - planeCenter_;
		if (meshCenterDir * planeNor_ > 0)
			planeNor_ = -planeNor_;
	}
	planeSideNor_ = (planeMidNor_ ^ planeNor_).normalized();
	node_->addDrawable(SmartDrawable(new LineGeometry(planeCenter_, planeCenter_ + planeNor_ * 10.0f, Color4f::Red, 10.0f)));
	node_->addDrawable(SmartDrawable(new LineGeometry(planeCenter_, planeCenter_ + planeMidNor_ * 10.0f, Color4f::Green, 10.0f)));
	node_->addDrawable(SmartDrawable(new LineGeometry(planeCenter_, planeCenter_ + planeSideNor_ * 10.0f, Color4f::Yellow, 10.0f)));

}

bool CropScanRod::PImpl::bIntersection(Vec3 rayOrigin, Vec3 rayNormal)
{
	ray_.SetOrigin(Vec3d::Construct(rayOrigin));
	ray_.SetDirection(Vec3d::Construct(rayNormal));
	if (mq_->intersects(ray_, out_intersection_))
	{
		return true;
	}
	else
	{
		ray_.SetDirection(Vec3d::Construct(-rayNormal));
		if (mq_->intersects(ray_, out_intersection_))
			return true;
		else
			return false;
	}
}


std::set<TVertex*> CropScanRod::PImpl::findSeedPts()
{
#if 0
	std::set<TVertex*> seedVertex;
	bool isSeed;
	for (auto& iter : pMesh_->face)
	{
		if((iter.IsB(0)|| iter.IsB(1)|| iter.IsB(2))&&)
			continue;
		//if (acos(iter.N() * planeNor_)<M_PI/4)
		{
			std::unordered_set<TFace*> ringfaces, traverseFaces;
			ringfaces.insert(&iter);
			for (int i = 0; i < 2; i++)
			{
				auto faceIter = ringfaces.begin();
				int range = ringfaces.size();
				for (int j = 0; j < range; j++)
				{
					if (traverseFaces.find(*faceIter) != traverseFaces.end())
					{
						++faceIter;
						continue;
					}
					else
					{
						ringfaces.insert((*faceIter)->cFFp(0));
						ringfaces.insert((*faceIter)->cFFp(1));
						ringfaces.insert((*faceIter)->cFFp(2));

						traverseFaces.insert(*faceIter);
					}
					++faceIter;
				}
			}
			isSeed = true; 
			for (auto& it : ringfaces)
			{
				float nor_angle = acos(it->N() * iter.N());
				if (nor_angle > M_PI/72)
				{
					isSeed = false;
					continue;
				}
			}
			if (isSeed)
			{
				seedVertex.insert(iter.V(0));
				seedVertex.insert(iter.V(1));
				seedVertex.insert(iter.V(2));
			}
		}
	}

	//筛选种子点
	Vec3Vector seedFaces;
	Vec3Vector ringVertexs;
	IndexDistVector out_result;
	for (auto iter : seedVertex)
	{
		if(iter->IsB())
			continue;
		if (bIntersection(iter->P(), iter->N()))
			continue;
		ann_->queryRadius(iter->P(), 2.2f, out_result);
		FaceList ringFaces;
		for (auto& itor : out_result)
		{
			TVertex* queryP = pMesh_->getVertex(itor.first);
			float theta = acos(queryP->N() * iter->N());
			if (theta > M_PI * 0.48 && theta < M_PI * 0.52)
			{
				isSeed = true;
				pMesh_->getRingFace(queryP, ringFaces);
				for (auto it : ringFaces)
				{
					float nor_angle = acos(it->N() * queryP->N());
					if (nor_angle > M_PI / 36)
					{
						isSeed = false;
						continue;
					}
				}
				if (isSeed)
				{
					ringVertexs.push_back(queryP->P());
					seedFaces.push_back(iter->P());
				}
				else
					continue;
			}
		}
	}
#endif
	
	std::map<TVertex*, float> pointProjectDist;
	Vec3 projectedPoint;
	float proDist;
	Vec3 meshCenter = pMesh_->getBoundingBox().Center();
	for (auto& iter:pMesh_->vert)
	{
		GeomOp::calcProjectPointOnLine(iter.P(),meshCenter, meshCenter + planeNor_ * 10.f, projectedPoint);
		if ((projectedPoint - meshCenter) * planeNor_ > 0)
		{
			proDist = (projectedPoint - meshCenter).Length();
		}
		else
		{
			proDist = -(projectedPoint - meshCenter).Length();
		}
		//pointProjectDist[&iter] = proDist;
		pointProjectDist.insert(std::make_pair(&iter, proDist));
	}

	//对值排序
	std::vector<std::pair<TVertex*, float>> sortProjectDist;
	for (auto it = pointProjectDist.begin(); it != pointProjectDist.end(); it++) {
		sortProjectDist.push_back(std::pair<TVertex*, float>(it->first, it->second));
	}
	sort(sortProjectDist.begin(), sortProjectDist.end(), orderValue);


	SelectFilter::clearSelect(pMesh_.get());
	Vec3Vector radiusQueryList;//show-----------
	std::map<TVertex*, std::vector<TVertex*>> perSeedQueryList;
	std::map<TVertex*, BoundingBox> queryPtsBox;
	
	std::set<TVertex*> seedPts;
	std::set<int> seedRadiusQueryPts;
	IndexDistVector queryPts;
	auto iter = --sortProjectDist.end();
	seedPts.insert(iter->first);
	iter->first->SetV();
	ann_->queryRadius(iter->first->P(), queryRadius_, queryPts);
	for (auto& it : queryPts)
	{
		auto cutTVertex = pMesh_->getVertex(it.first);
		cutTVertex->SetV();
		perSeedQueryList[iter->first].push_back(cutTVertex);
		radiusQueryList.push_back(cutTVertex->P());
	}
	for (; iter != sortProjectDist.begin();--iter)
	{
		if(iter->first->IsV()/* || iter->first->IsB()*/)
			continue;
		bool bSeed_ = false;//种植点的距离是否太近
		for (auto seed : seedPts)
		{
			float seed_length = (iter->first->P() - seed->P()).Length();
			if (seed_length < queryRadius_*2)
			{
				bSeed_ = true;
				break;
			}
		}
		if (bSeed_)
		{
			continue;
		}
		iter->first->SetV();
		//radiusQueryList.push_back(iter->first->P());
		ann_->queryRadius(iter->first->P(), queryRadius_, queryPts);
		seedPts.insert(iter->first);
		for (auto& it : queryPts)
		{
			auto cutTVertex = pMesh_->getVertex(it.first);
			cutTVertex->SetV();
			perSeedQueryList[iter->first].push_back(cutTVertex);
			radiusQueryList.push_back(cutTVertex->P());
		}
		SmartMesh copy_mesh = Mesh::copyMesh(*pMesh_.get());
		if (seedPts.size() >= scanRodNum_)
		{
			std::map<TVertex*, shared_ptr<ObbFilter>> perSeedQueryBox;
			for (auto ps : perSeedQueryList)
			{
				SelectFilter::clearSelect(copy_mesh.get());
				for (auto sq : ps.second)
				{
					copy_mesh->getVertex(sq->id())->SetS();
				}
				SelectFilter::selectFacesFromVertices(copy_mesh.get(), false);
				SmartMesh per_mesh(new Mesh);
				Mesh::appendMesh(*per_mesh, *copy_mesh,true);
				obbObject_ = std::make_shared<ObbFilter>(per_mesh.get());
				
				perSeedQueryBox[ps.first] = obbObject_;
			}
			for (auto pbox : perSeedQueryBox)
			{
				float box_d = pbox.second->depth();
				float box_w = pbox.second->width();
				float box_h = pbox.second->height();
				if (pbox.second->depth() < queryRadius_*0.8|| pbox.second->width() < queryRadius_ * 0.8|| pbox.second->height() < queryRadius_ * 0.8)
				{
					auto it1 = perSeedQueryBox.find(pbox.first);
					perSeedQueryBox.erase(it1);
					perSeedQueryList.erase(it1->first);
					//auto it2 = seedPts.find(pbox.first);
					seedPts.erase(pbox.first);
				}
			}
			if(seedPts.size() >= scanRodNum_)
				break;
		}
	}

	pointProjectDist.clear();
	for (auto ss : sortProjectDist)
	{
		pointProjectDist[ss.first] = ss.second;
	}

	//种子点更新到planeNor_方向上的最高点
	std::set<TVertex*> newSeedPts;
	for (auto iter : perSeedQueryList)
	{
		TVertex* newSeed = iter.first;
		for (auto itor : iter.second)
		{
			if (pointProjectDist[iter.first] < pointProjectDist[itor])
			{
				newSeed = itor;
			}
		}
		newSeedPts.insert(newSeed);
	}

	Mesh::VertexAttributeByte vertexAttr;
	pMesh_->createVertexAttribute("isSeed", vertexAttr);
	pMesh_->setVertexAttributeValue(vertexAttr, 0);
	for (auto itor : newSeedPts)
	{
		vertexAttr[itor] = 1;
	}

	//node_->addDrawable(pMesh_);
	Vec3Vector seedList;
	for (auto itor : newSeedPts)
	{
		seedList.push_back(itor->P());
	}
	
	//node_->addDrawable(SmartDrawable(new PointGeometry(radiusQueryList, Color4f::Blue, 5.0f)));
	node_->addDrawable(SmartDrawable(new PointGeometry(seedList, Color4f::Red, 10.0f)));
	//node_->addDrawable(pMesh_);
	return newSeedPts;
}

void CropScanRod::PImpl::tracerseCrop2(std::set<TVertex*>& seedPts)
{
	std::map<TVertex*, ProPts> vPerDirProPoint;//每个点在三个方向上的投影位置
	Vec3 sideProPoint, midProPoint, norProPoint;
	Vec3 meshCenter = pMesh_->getBoundingBox().Center();
	for (auto& iter : pMesh_->vert)
	{
		GeomOp::calcProjectPointOnLine(iter.P(), meshCenter, meshCenter + planeNor_ * 100.f, norProPoint);
		GeomOp::calcProjectPointOnLine(iter.P(), meshCenter, meshCenter + planeSideNor_ * 100.f, sideProPoint);
		GeomOp::calcProjectPointOnLine(iter.P(), meshCenter, meshCenter + planeMidNor_ * 100.f, midProPoint);

		vPerDirProPoint[&iter].midP = midProPoint;
		vPerDirProPoint[&iter].norP = norProPoint;
		vPerDirProPoint[&iter].sideP = sideProPoint;
	}

	std::vector<SmartMesh> rodMeshList;
	vector<float> curvatureList;
	SelectFilter::clearSelect(pMesh_.get());
	std::set<TVertex*> curSeedTraverseList;
	std::map<TVertex*, std::set<TVertex*>> allSeedTraverseList;
	Vec3Vector trarvalList;//show
	Vec3Vector visionList;//show
	std::set<TVertex*> newSeedPts;
	for (auto seed : seedPts)
	{
		for (int i = 0; i < pMesh_->numVertex(); i++)
		{
			pMesh_->getVertex(i)->ClearV();
		}
		SelectFilter::clearSelect(pMesh_.get());
		curSeedTraverseList.clear();
		std::queue<TVertex*> selectP;
		selectP.push(seed);
		seed->SetV();
		seed->SetS();
		visionList.push_back(seed->P());
		while (!selectP.empty())
		{
			auto it = selectP.front();
			TVertex* curTVertex = it;
			selectP.pop();
			VertexList ringVertexList;
			pMesh_->getRingVertex(curTVertex, ringVertexList);
			for (auto ringV : ringVertexList)
			{
				if (ringV->IsV())
					continue;
				float seedDist = (ringV->P() - seed->P()).Length();
				float midDist = (vPerDirProPoint[ringV].midP - vPerDirProPoint[seed].midP).Length();
				float norDist = (vPerDirProPoint[ringV].norP - vPerDirProPoint[seed].norP).Length();
				float sideDist = (vPerDirProPoint[ringV].sideP - vPerDirProPoint[seed].sideP).Length();
				if (seedDist>= scanCrodHeight_|| norDist >= scanCrodHeight_
					|| midDist>=queryRadius_|| sideDist >= queryRadius_)//
				{
					ringV->SetV();
					visionList.push_back(ringV->P());
					break;
				}
				else
				{
					ringV->SetV();
					ringV->SetS();
					selectP.push(ringV);
					visionList.push_back(ringV->P());
				}
			}
			curSeedTraverseList.insert(curTVertex);
		}
		allSeedTraverseList[seed] = curSeedTraverseList;

		SelectFilter::selectFacesFromVertices(pMesh_.get(), false);
		SmartMesh per_mesh(new Mesh);
		Mesh::VertexAttributeByte vertexAttr;
		per_mesh->createVertexAttribute("isSeed", vertexAttr);
		//per_mesh->setVertexAttributeValue(vertexAttr, 0);
		Mesh::appendMesh(*per_mesh, *pMesh_, true);
		rodMeshList.push_back(per_mesh);
	}

	int j = 0;
	for (auto iter : rodMeshList)
	{
		//node_->addDrawable(iter);
		std::stringstream ss;
		ss << "D:/ProjectTestData/CropScanRod/prepare/scan_rod_";
		ss << j;
		ss << ".um";
		MeshIO::writeMeshd(*iter, ss.str().c_str());
		j++;
	}



	Vec3Vector topList;//show--
	int l = 0;
	Vec3Vector travelList;
	std::vector<SmartMesh> resultMeshList;
	Vec3Vector showPts1,showPts2;
	for (auto iter : rodMeshList)
	{
		CurvatureFilter::MeanAndGaussian(*iter.get());
		Mesh::VertexAttributeByte vertexAttr;
		iter->getVertexAttribute("isSeed", vertexAttr);
		TVertex* seedVertex;
		for (auto& it:iter->vert)
		{
			if (vertexAttr[it] == 1)
			{
				seedVertex = &it;
				break;
			}
		}
		
		TVertex* featureVertex;
		SelectFilter::clearSelect(iter.get());
		std::queue<TVertex*> selectPts;
		selectPts.push(seedVertex);
		seedVertex->SetV();
		seedVertex->SetS();
		while (!selectPts.empty())
		{
			auto curVertex = selectPts.front();
			selectPts.pop();
			VertexList ringVertexList;
			iter->getRingVertex(curVertex, ringVertexList);
			bool loop_ = true;
			for (auto rv : ringVertexList)
			{
				if(rv->IsV())
					continue;
				rv->SetV();
				rv->SetS();
				if (rv->Kh() <0)
				{
					showPts1.push_back(rv->P());
				}
				else
				{
					showPts2.push_back(rv->P());
					selectPts.push(rv);
				}
			}
		}

		SelectFilter::selectFacesFromVertices(iter.get(), false);
		SmartMesh per_mesh(new Mesh);
		Mesh::VertexAttributeByte vertexAttr1;
		per_mesh->createVertexAttribute("isSeed", vertexAttr1);
		Mesh::appendMesh(*per_mesh, *iter, true);
		per_mesh->dirtyMesh();
		CleanFilter::removeSmallComponents(per_mesh.get());
		obbObject_ = std::make_shared<ObbFilter>(per_mesh.get());
		//node_->addDrawable(SmartDrawable(new LineGeometry(obbObject_->getCenter(), obbObject_->getCenter() + obbObject_->getAxes().main * 10.f, Color4f::Red, 10.0f)));
		//node_->addDrawable(SmartDrawable(new LineGeometry(obbObject_->getCenter(), obbObject_->getCenter() + obbObject_->getAxes().second * 10.f, Color4f::Green, 10.0f)));
		//node_->addDrawable(SmartDrawable(new LineGeometry(obbObject_->getCenter(), obbObject_->getCenter() + obbObject_->getAxes().third * 10.f, Color4f::Blue, 10.0f)));
		FillHole fh;
		fh.set_mesh(per_mesh);
		std::vector<IntVector> all_hole = fh.detect_hole();
		int max_size = -1, max_index = -1,hole_size = -1;
		for (int i =0;i<all_hole.size();i++)
		{
			hole_size = all_hole[i].size();
			if (hole_size > max_size)
			{
				max_size = hole_size;
				max_index = i;
			}
		}
		Vec3 boundaryCenter = Vec3{ 0.0f,0.0f,0.0f };
		for (int bid : all_hole[max_index])
		{
			boundaryCenter += (per_mesh->getVertex(bid)->P());
		}
		boundaryCenter /= all_hole[max_index].size();
		node_->addDrawable(SmartDrawable(new PointGeometry(obbObject_->getCenter(), Color4f::Yellow, 10.0f)));
		node_->addDrawable(SmartDrawable(new PointGeometry(boundaryCenter, Color4f::Yellow, 10.0f)));

		Vec3 center_dir = (obbObject_->getCenter() - boundaryCenter).normalized();

		VertexList topPts;//show--
		for (int i = 0; i < iter->numVertex(); i++)
		{
			iter->getVertex(i)->ClearV();
		}
		SelectFilter::clearSelect(iter.get());
		selectPts.push(seedVertex);
		seedVertex->SetV();
		seedVertex->SetS();
		while (!selectPts.empty())
		{
			auto curVertex = selectPts.front();
			selectPts.pop();
			VertexList ringVertexList;
			iter->getRingVertex(curVertex, ringVertexList);
			for (auto rv : ringVertexList)
			{
				if (rv->IsV())
					continue;
				if (acos(rv->N()* center_dir)>M_PI/4 /*&& (rv->P()-obbObject_->getCenter()) * center_dir>0*/)
				{
					rv->SetV();
					//break;
				}
				else
				{
					rv->SetV();
					rv->SetS();
					topPts.push_back(rv);
					selectPts.push(rv);
				}
			}
		}
		// 計算頂部中心
		Vec3 topCenter = Vec3{ 0.0f,0.0f,0.0f };
		Vec3 topNor = Vec3{ 0.0f,0.0f,0.0f };
		if (topPts.empty())
		{
			GeomOp::calcProjectPointOnLine(seedVertex->P(), obbObject_->getCenter(), boundaryCenter, topCenter);
			topNor = (topCenter - boundaryCenter).normalized();
		}
		else
		{
			for (auto tv : topPts)
			{
				topList.push_back(tv->P());
				topCenter += tv->P();
				topNor += tv->N();
			}
			topNor = topNor.normalized();
			topCenter /= topPts.size();
		}
		
		node_->addDrawable(SmartDrawable(new LineGeometry(topCenter, topCenter + topNor*20.f, Color4f::Red, 10.0f)));
		node_->addDrawable(SmartDrawable(new PointGeometry(topCenter, Color4f::Green, 10.0f)));
 
		VertexList boundaryVertexList;
		iter->getBoundaryVertexList(boundaryVertexList);
		for (int i = 0; i < iter->numVertex(); i++)
		{
			iter->getVertex(i)->ClearV();
		}
		SelectFilter::clearSelect(iter.get());
		VertexList deleteVertex;
		for (auto bv : boundaryVertexList)
		{
			/*if(bv->N() * topNor < 0)
				continue;*/
			selectPts.push(bv);
			deleteVertex.push_back(bv);
			bv->SetV();
		}
		while (!selectPts.empty())
		{
			auto curVertex = selectPts.front();
			selectPts.pop();
			VertexList ringVertexList;
			iter->getRingVertex(curVertex, ringVertexList);
			for (auto rv : ringVertexList)
			{
				if (rv->IsV())
					continue;
				rv->SetV();
				if(acos(rv->N() * topNor)<M_PI*0.4&&(rv->P()-obbObject_->getCenter())*topNor<0)
				{
					//rv->SetV();
					//rv->SetS();
					deleteVertex.push_back(rv);
					travelList.push_back(rv->P());
					selectPts.push(rv);
				}
			}
		}
		for (auto dv:deleteVertex)
		{
			iter->deleteVertex(dv);
		}
		iter->compact();
		iter->dirtyMesh();
		CleanFilter::removeSmallComponents(iter.get());
		resultMeshList.push_back(iter);

#if 0
		selectPts.push(seedVertex);
		seedVertex->SetV();
		seedVertex->SetS();
		while (!selectPts.empty())
		{
			auto curVertex = selectPts.front();
			selectPts.pop();
			VertexList ringVertexList;
			iter->getRingVertex(curVertex, ringVertexList);
			bool loop_ = true;
			for (auto rv : ringVertexList)
			{
				if (rv->IsV())
					continue;
				float pro_dist = GeomOp::calcDistFromPointToLine(rv->P(), boundaryCenter, obbObject_->getCenter());
				if (pro_dist > queryRadius_*0.75)
				{
					rv->SetV();
				}
				else
				{
					rv->SetV();
					rv->SetS();
				}
				selectPts.push(rv);
			}
		}
		SelectFilter::selectFacesFromVertices(iter.get(), false);
		SmartMesh per_mesh2(new Mesh);
		/*per_mesh2->createVertexAttribute("isSeed", vertexAttr);
		per_mesh2->setVertexAttributeValue(vertexAttr, 0);*/
		Mesh::appendMesh(*per_mesh2, *iter, true);
		per_mesh2->dirtyMesh();
		//CleanFilter::removeSmallComponents(per_mesh2.get());
		resultMeshList.push_back(per_mesh2);
#endif
	}
	node_->addDrawable(SmartDrawable(new PointGeometry(topList, Color4f::Green, 1.0f)));
	node_->addDrawable(SmartDrawable(new PointGeometry(travelList, Color4f::Blue, 1.0f)));

	int k = 0;
	for (auto iter : resultMeshList)
	{
		node_->addDrawable(iter);
		std::stringstream ss;
		ss << "D:/ProjectTestData/CropScanRod/traveralMesh/scan_rod_";
		ss << k;
		ss << ".um";
		MeshIO::writeMeshd(*iter, ss.str().c_str());
		k++;
	}

	//node_->addDrawable(SmartDrawable(new PointGeometry(showPts1, Color4f::Red, 3.0f)));
	//node_->addDrawable(SmartDrawable(new PointGeometry(showPts2, Color4f::Green, 3.0f)));

#if 0 
	Vec3Vector showPts;
	Vec3Vector showCurvature;
	VertexList deleteVertexList;
	VertexList boundaryVertexList;
	int i = 0;
	for (auto iter : rodMeshList)
	{
		iter->dirtyMesh();
		//iter->enableVertexCurvature(true);
		//iter->enableVertexCurvatureDir(true);
		CurvatureFilter::MeanAndGaussian(*iter.get());
		/*float minCurvature = FLT_MAX;
		float maxCurvature = FLT_MIN;
		for (auto& vc : iter->vert)
		{
			if (vc.Kh() < minCurvature)
			{
				minCurvature = vc.Kh();
			}
			else if (vc.Kh() > maxCurvature)
			{
				maxCurvature = vc.Kh();
			}
		}
		float rangeLength = maxCurvature - minCurvature;*/

		/*for (auto& vc : iter->vert)
		{
			float curCurvature = (vc.Kh() - minCurvature) / rangeLength;
			if (curCurvature > 0.8 || vc.Kh() <0 )
				showCurvature.push_back(vc.P());
		}*/

		boundaryVertexList.clear();
		deleteVertexList.clear();
		iter->getBoundaryVertexList(boundaryVertexList);
		std::queue<TVertex*> startPts;
		for (auto bv : boundaryVertexList)
		{
			startPts.push(bv);
		}
		while (!startPts.empty())
		{
			auto curVertex = startPts.front();
			startPts.pop();
			curVertex->SetV();
			VertexList ringVertex;
			iter->getRingVertex(curVertex, ringVertex);
			bool bDelete_ = false;
			for (auto rv : ringVertex)
			{
				if (rv->IsV())
					continue;
				if (rv->Kh() > 0.2/*rv->Kh()>0||rv->Kg()<0*/)
				{
					//bDelete_ = true;
					deleteVertexList.push_back(curVertex);
					deleteVertexList.push_back(rv);
					//rv->SetV();
					//startPts.push(rv);
					continue;
					//break;
				}
				rv->SetV();
				startPts.push(rv);
			}
		}

		

		//删除面片
		for (auto dv : deleteVertexList)
		{
			iter->deleteVertex(dv);
		}
		iter->compact();
		iter->dirtyMesh();
		std::cout << iter->numFace() << endl;
		std::cout << deleteVertexList.size() << endl;
		CleanFilter::removeSmallComponents(iter.get());




		//扫描杆中心
		/*Vec3 seedCenter = Vec3{ 0.f,0.f,0.f };
		int innnerNum = 0;
		mq_ = iter->createTree(ITree::TT_QUERY)->asQuery();
		for (auto& vc : iter->vert)
		{
			if (vc.IsV())
				continue;
			vc.SetV();
			VertexList ringVertex;
			iter->getRingVertex(&vc, ringVertex);
			Vec3 dir = Vec3{0.f,0.f,0.f};
			bool bJump_ = false;
			for (auto rv : ringVertex)
			{
				if (rv->IsV())
				{
					bJump_ = true;
					break;
				}
				dir += rv->N();
				rv->SetV();
			}
			if (bJump_)
				continue;
			ray_.SetOrigin(Vec3d::Construct(vc.P()));
			ray_.SetDirection(Vec3d::Construct(vc.N()));
			if (mq_->intersects(ray_, out_intersection_))
			{
				showPts.push_back(vc.P());
				seedCenter += vc.P();
				innnerNum++;
			}
		}
		seedCenter /= innnerNum;
		obbObject_ = std::make_shared<ObbFilter>(iter.get());
		float d = abs(obbObject_->depth()- queryRadius_ / 2);
		float w = abs(obbObject_->width() - queryRadius_ / 2);
		float h = abs(obbObject_->height() - queryRadius_ / 2);
		Vec3 seedDir = h > (d > w ? d : w) ? obbObject_->getAxes().second : (d > w ? obbObject_->getAxes().third : obbObject_->getAxes().main);

		node_->addDrawable(SmartDrawable(new LineGeometry(seedCenter, seedCenter + seedDir *10.f, Color4f::Green, 2.0f)));*/


		iter->setBack(true);
		node_->addDrawable(iter);
		//node_->addDrawable(SmartDrawable(new PointGeometry(showCurvature, Color4f::Red, 1.0f)));
		std::stringstream ss;
		ss << "D:/ProjectTestData/CropScanRod/output/scan_rod_";
		ss << i;
		ss << ".um";
		MeshIO::writeMeshd(*iter, ss.str().c_str());
		i++;

	}
	//node_->addDrawable(pMesh_);
	node_->addDrawable(SmartDrawable(new PointGeometry(showPts, Color4f::Red, 5.0f)));
#endif

#if 0
	FaceList deleteFace;
	FaceList boundaryFaceList;
	int i = 0;
	for (auto iter : rodMeshList)
	{
		boundaryFaceList.clear();
		deleteFace.clear();
		for (auto it : iter->face)
		{
			if (it.IsB(0)|| it.IsB(1)|| it.IsB(2))
			{
				boundaryFaceList.push_back(&it);
			}
		}

		std::queue<TFace*> startPts;
		for (auto it : boundaryFaceList)
		{
			startPts.push(it);
			it->SetV();
		}
		while (!startPts.empty())
		{
			auto curFace = startPts.front();
			startPts.pop();
			FaceList ringFace;
			iter->getRingFace(curFace, ringFace);
			bool bDelete_ = false;
			for (auto rf : ringFace)
			{
				if(rf->IsV())
					continue;
				float angle = acos(rf->N() * curFace->N());
				if (angle > M_PI/16)
				{
					bDelete_ = true;
					startPts.push(rf);
					continue;
				}
			}
			if (bDelete_)
			{
				deleteFace.push_back(curFace);
				curFace->SetV();
			}
			else
			{
				curFace->SetV();
			}
		}

		//删除面片
		for (auto it : deleteFace)
		{
			iter->deleteFaceAndVertex(it);
		}
		cout << iter->numFace() << endl;
		cout << deleteFace.size() << endl;
		iter->compact();
		iter->dirtyMesh();
		node_->addDrawable(iter);

		std::stringstream ss;
		ss << "D:/ProjectTestData/CropScanRod/output/scan_rod_";
		ss << i;
		ss << ".um";
		MeshIO::writeMeshd(*iter, ss.str().c_str());
		i++;
	}
}
#endif
}
void CropScanRod::PImpl::tracerseCrop(std::set<TVertex*> seedPts)
{
	//计算曲率
	CurvatureFilter::MeanAndGaussian(*pMesh_.get());
	vector<float> curvatureList;
	SelectFilter::clearSelect(pMesh_.get());
	std::set<TVertex*> curSeedQueryList;
	std::map<TVertex* ,std::set<TVertex*>> allSeedQueryList;
	Vec3Vector trarvalList;
	Vec3Vector visionList;//show
	for (auto seed : seedPts)
	{
		if (seed->IsV())
			continue;
		SelectFilter::clearSelect(pMesh_.get());
		curSeedQueryList.clear();
		std::queue<TVertex*> selectP;
		selectP.push(seed);
		seed->SetV();
		seed->SetS();
		visionList.push_back(seed->P());
		while (!selectP.empty())
		{
			auto it = selectP.front();
			TVertex* curTVertex = it;
			selectP.pop();
			VertexList ringVertexList;
			pMesh_->getRingVertex(curTVertex, ringVertexList);
			for (auto ringV : ringVertexList)
			{
				if (ringV->IsV())
					continue;
				float seedDist = (ringV->P() - seed->P()).Length();
				/*if (ringV->Kh() < -1)*/
				if(seedDist >=scanCrodHeight_)
				{
					ringV->SetV();
					visionList.push_back(ringV->P());
					break;
				}
				else
				{
					ringV->SetV();
					ringV->SetS();
					selectP.push(ringV);
					visionList.push_back(ringV->P());
				}
			}
			curSeedQueryList.insert(curTVertex);
		}
		for (auto iter : curSeedQueryList)
		{
			trarvalList.push_back(iter->P());
		}
		allSeedQueryList[seed] = curSeedQueryList;
	}
	node_->addDrawable(SmartDrawable(new PointGeometry(visionList, Color4f::Blue, 5.0f)));
	//node_->addDrawable(SmartDrawable(new PointGeometry(trarvalList, Color4f::Green, 5.0f)));
	//归一化每一个遍历点集的曲率
	//std::vector<std::vector<float>> allSeedCurvatureList;
	//std::vector<float> perSeedCurvatureList;
	trarvalList.clear();
	for (int  i =0;i<pMesh_->numVertex();i++)
	{
		pMesh_->getVertex(i)->ClearV();
	}
	for (auto iter : allSeedQueryList)
	{
		auto seed_ = iter.first;
		if (seed_->IsV())
			continue;
		SelectFilter::clearSelect(pMesh_.get());
		float minCurvature = FLT_MAX;
		float maxCurvature = FLT_MIN;
		for (auto it : iter.second)
		{
			if (it->Kh() < minCurvature)
			{
				minCurvature = it->Kh();
			}
			else if(it->Kh()> maxCurvature)
			{
				maxCurvature = it->Kh();
			}
		}
		float rangeLength = maxCurvature - minCurvature;

		//重新遍历裁剪扫描杆
		
		std::queue<TVertex*> selectP;
		curSeedQueryList.clear();
		seed_->SetV();
		seed_->SetS();
		selectP.push(seed_);
		while (!selectP.empty())
		{
			auto it = selectP.front();
			TVertex* curTVertex = it;
			selectP.pop();
			VertexList ringVertexList;
			pMesh_->getRingVertex(curTVertex, ringVertexList);
			for (auto ringV : ringVertexList)
			{
				if (ringV->IsV()|| iter.second.find(ringV) == iter.second.end())
					continue;
				float seedDist;
				if ((ringV->P() - seed_->P()) * seed_->N() < 0)
				{
					seedDist = abs((ringV->P() - seed_->P()) * seed_->N());
				}
				else
				{
					seedDist = (ringV->P() - seed_->P()).Length();
				}
				//float cone_angle = acos(ringV->N() * curTVertex->N());
				//cout << "cone_angle------" << cone_angle << endl;
				//if (cone_angle > 1)
				//if(seedDist >= scanCrodHeight_)
				float curCurvature = (ringV->Kh() - minCurvature) / rangeLength;
				//cout << "curCurvature------" << curCurvature << endl;
				if (curCurvature < 0.1)
				{
					ringV->SetV();
					visionList.push_back(ringV->P());
					break;
				}
				else
				{
					ringV->SetV();
					ringV->SetS();
					selectP.push(ringV);
					visionList.push_back(ringV->P());
				}
			}
			curSeedQueryList.insert(curTVertex);
		}
		for (auto iter : curSeedQueryList)
		{
			trarvalList.push_back(iter->P());
		}
		allSeedQueryList[seed_] = curSeedQueryList;
	}
	node_->addDrawable(SmartDrawable(new PointGeometry(trarvalList, Color4f::Green, 5.0f)));
	
}

bool CropScanRod::PImpl::orderValue(std::pair<TVertex*, float> a, std::pair<TVertex*, float> b)
{
	return a.second < b.second;
}

bool CropScanRod::PImpl::apply()
{
	calcOcclusalSurface();
	std::set<TVertex*> seedPts = findSeedPts();
	tracerseCrop2(seedPts);
	
	//SurfaceMeshHoleFiller::fillSomeHoles(pMesh_.get(), 0, 0, 2, true);
	//MeshIO::writeMesh(*pMesh_, "CropScanRod/ScanRod.off");
	//std::vector<double> sdf_list1, sdf_list2, sdf_list3;
	// // create and read Polyhedron
	//SM mesh;
	//if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(CGAL::data_file_path("CropScanRod/ScanRod.off"), mesh) ||
	//	!CGAL::is_triangle_mesh(mesh))
	//{
	//	std::cerr << "Invalid input file." << std::endl;
	//	return EXIT_FAILURE;
	//}

	//typedef SM::Property_map<face_descriptor, double> Facet_double_map;
	//Facet_double_map sdf_property_map;
	//sdf_property_map = mesh.add_property_map<face_descriptor, double>("f:sdf").first;

	//auto start_time_ = std::chrono::system_clock::now();
	//CGAL::sdf_values(mesh, sdf_property_map);
	//typedef SM::Property_map<face_descriptor, std::size_t> Facet_int_map;
	//Facet_int_map segment_property_map = mesh.add_property_map<face_descriptor, std::size_t>("f:sid").first;
	//std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);
	//const std::size_t number_of_clusters = 15;
	//const double smoothing_lambda = 0.3;
	//CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);
	//auto end_time = std::chrono::system_clock::now();
	//double totalTime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time_).count() / 1e3;
	//std::cout << "[info] all time is: " << totalTime << " ms\n";
	//
	//typedef CGAL::Face_filtered_graph<SM> Filtered_graph;
	////print area of each segment and then put it in a Mesh and print it in an OFF file
	//Filtered_graph segment_mesh(mesh);
	//for (std::size_t id = 0; id < number_of_segments; ++id)
	//{
	//	segment_mesh.set_selected_faces(id, segment_property_map);
	//	if(segment_mesh.number_of_faces()<100)
	//		continue;
	//	std::cout << "Segment " << id << "'s area is : " << CGAL::Polygon_mesh_processing::area(segment_mesh) << std::endl;
	//	SM out;
	//	CGAL::copy_face_graph(segment_mesh, out);
	//	std::ostringstream oss;
	//	oss << "CropScanRod/Segment_" << id << ".off";
	//	std::ofstream os(oss.str().data());
	//	os << out;
	//}

	//for (face_descriptor f : faces(mesh))
	//{
	//	sdf_list1.push_back(sdf_property_map[f]);
	//	//std::cout << sdf_property_map[f] << " ";
	//}
	////std::cout << std::endl;

	/*node_->addDrawable(pMesh_);
	pMesh_->setColorMode(CMPerVert);


	std::size_t i = 0;
	for (Polyhedron::Facet_iterator it = mesh.facets_begin(); it != mesh.facets_end(); ++it, ++i)
	{
		it->id() = i;
	}

	for (face_descriptor f : faces(mesh))
	{
		double cur_sdf = sdf_property_map[f];
		auto idx = f->id();
		TFace* curFace = pMesh_->getFace(f->id());
		if (cur_sdf < 0.25)
		{
			curFace->V(0)->C() = Color4b(255, 0, 0, 0);
			curFace->V(1)->C() = Color4b(255, 0, 0, 0);
			curFace->V(2)->C() = Color4b(255, 0, 0, 0);
		}
		else if(cur_sdf < 0.5)
		{
			curFace->V(0)->C() = Color4b(0, 255, 0, 0);
			curFace->V(1)->C() = Color4b(0, 255, 0, 0);
			curFace->V(2)->C() = Color4b(0, 255, 0, 0);
		}
		else if (cur_sdf < 0.75)
		{
			curFace->V(0)->C() = Color4b(0, 0, 255, 0);
			curFace->V(1)->C() = Color4b(0, 0, 255, 0);
			curFace->V(2)->C() = Color4b(0, 0, 255, 0);
		}
		else
		{
			curFace->V(0)->C() = Color4b(0, 255, 255, 0);
			curFace->V(1)->C() = Color4b(0, 255, 255, 0);
			curFace->V(2)->C() = Color4b(0, 255, 255, 0);
		}
	}*/
	//// ------------------------------------------------------------------------------------
	//double boxRadius = 0.05;
	//Vec3d face_centroid;
	//SDF_property_map face_property;
	//MeshQuery::IntersectionList out_intersection_list;
	//sdf_map.clear();
	//mq_ = pMesh_->createTree(ITree::TT_QUERY)->asQuery();
	//for (auto& iter : pMesh_->face)
	//{
	//	face_centroid = Vec3d::Construct((iter.V(0)->P() + iter.V(1)->P() + iter.V(2)->P()) / 3);
	//	ray_.SetOrigin(face_centroid);
	//	ray_.SetDirection(Vec3d::Construct(-iter.N()).normalized());
	//	if (mq_->intersects(ray_, out_intersection_list))
	//	{
	//		/*Vec3d interP1 = out_intersection_list[0].intersection_point;
	//		ray_.SetDirection(Vec3d::Construct(-iter.N()).normalized());
	//		if (mq_->intersects(ray_, out_intersection_list))
	//		{
	//			Vec3d interP2 = out_intersection_list[0].intersection_point;
	//			face_property.sdf_value = (interP1 - interP2).Length() - boxRadius * 2;
	//			face_property.postprocess = true;
	//			if(face_property.sdf_value<0)
	//				face_property.postprocess = false;
	//		}
	//		else
	//		{
	//			face_property.postprocess = false;
	//		}*/
	//		Vec3d interP = out_intersection_list[0].intersection_point;
	//		Vec3 interPlaneNor = Vec3::Construct(pMesh_->getFace(out_intersection_list[0].triangle_index)->N()).normalized();
	//		float theta = asin((-iter.N().normalized() ^ interPlaneNor).Length());
	//		face_property.sdf_value = (interP - face_centroid).Length() -(1 / cos(theta) - 1) * boxRadius;
	//		face_property.postprocess = true;
	//		if (face_property.sdf_value < 0)
	//			face_property.postprocess = false;
	//		face_property.postprocess = true;
	//	}
	//	else
	//	{
	//		face_property.postprocess = false;
	//	}
	//	sdf_map[&iter] = face_property;
	//}

	//FaceList ringFaceList;
	//std::vector<double> ringSdfList;
	//bool isMinSdf = false;
	//float avg_sdf_value = 0.0f,new_avg_sdf_value = 0.0f;
	//for (auto& iter : sdf_map)
	//{
	//	if(iter.second.postprocess)
	//		continue;
	//	pMesh_->getRingFace(iter.first, ringFaceList);
	//	for (auto& it : ringFaceList)
	//	{
	//		if (!sdf_map[it].postprocess)
	//		{
	//			isMinSdf = true;
	//			continue;
	//		}
	//		ringSdfList.push_back(sdf_map[it].sdf_value);
	//		avg_sdf_value += sdf_map[it].sdf_value;
	//	}
	//	std::sort(ringSdfList.begin(), ringSdfList.end());
	//	new_avg_sdf_value = avg_sdf_value;
	//	avg_sdf_value /= ringSdfList.size();
	//	if(isMinSdf)
	//	{
	//		for (auto& it : ringFaceList)
	//		{
	//			if (!sdf_map[it].postprocess)
	//			{
	//				new_avg_sdf_value += ringSdfList[0];
	//				sdf_map[it].sdf_value = ringSdfList[0];
	//				sdf_map[it].postprocess = true;
	//			}
	//		}
	//		iter.second.sdf_value = new_avg_sdf_value / ringFaceList.size();
	//		iter.second.postprocess = true;
	//	}
	//	else
	//	{
	//		iter.second.sdf_value = avg_sdf_value;
	//		iter.second.postprocess = true;
	//	}
	//}
	//BilateralSdfFilter(pMesh_);
	//for (auto& it : sdf_map)
	//{
	//	sdf_list2.push_back(it.second.sdf_value);
	//}

	//// ------------------------------------------------------------------------------------
	//sdf_map.clear();
	//avg_sdf_value = 0.0f;
	//new_avg_sdf_value = 0.0f;
	//SDF_property_map cur_sdf;
	//std::vector<double> ray_length;
	//mq_ = pMesh_->createTree(ITree::TT_QUERY)->asQuery();
	//for (auto& iter : pMesh_->face)
	//{
	//	face_centroid = Vec3d::Construct((iter.V(0)->P() + iter.V(1)->P() + iter.V(2)->P()) / 3);
	//	ray_.SetOrigin(face_centroid);
	//	MeshQuery::IntersectionList inter_;
	//	float theta = 2 * M_PI / numbers_of_rays_;
	//	Vec3 zDir = -iter.N(), xDir = { 1.0f,0,0 }, yDir = { 0,1.0f,0 };
	//	XOYOZ(zDir, xDir, yDir);
	//	for (int i = 0; i < numbers_of_rays_; i++)
	//	{
	//		float angle = i * theta;
	//		ray_.SetDirection(Vec3d::Construct((xDir * sqrt(3) * std::rand() * sin(angle) + yDir * sqrt(3) * std::rand() * cos(angle)).normalized()));
	//		if (mq_->intersects(ray_, inter_))
	//		{
	//			ray_length.push_back(inter_[0].dist);
	//		}
	//	}
	//	if (ray_length.size() == 0)
	//	{
	//		cur_sdf.postprocess = false;
	//	}
	//	else
	//	{
	//		double avg_sdk = 0.0;
	//		for (auto iter : ray_length)
	//		{
	//			avg_sdk += iter;
	//		}
	//		if (avg_sdk / ray_length.size()>0)
	//		{
	//			cur_sdf.postprocess = true;
	//			cur_sdf.sdf_value = avg_sdk / ray_length.size();
	//		}
	//		else
	//			cur_sdf.postprocess = false;
	//	}
	//	sdf_map[&iter] = cur_sdf;
	//}
	//for (auto& itor : sdf_map)
	//{
	//	if (itor.second.postprocess)
	//		continue;
	//	std::vector<TFace*> faceList;
	//	//遍历邻接面片
	//	pMesh_->getRingFace(itor.first, faceList);
	//	bool isMinValue = false;
	//	for (auto iter : faceList)
	//	{
	//		if (!sdf_map[iter].postprocess)
	//		{
	//			isMinValue = true;
	//			continue;
	//		}
	//		ringSdfList.push_back(sdf_map[iter].sdf_value);
	//		avg_sdf_value += sdf_map[iter].sdf_value;
	//	}
	//	new_avg_sdf_value = avg_sdf_value;
	//	avg_sdf_value /= ringSdfList.size();
	//	std::sort(ringSdfList.begin(), ringSdfList.end());
	//	if (isMinValue)
	//	{
	//		for (auto iter : faceList)
	//		{
	//			if (!sdf_map[iter].postprocess)
	//			{
	//				sdf_map[iter].sdf_value = ringSdfList[0];
	//				sdf_map[iter].postprocess = true;
	//				new_avg_sdf_value += ringSdfList[0];
	//			}
	//		}
	//		itor.second.sdf_value = new_avg_sdf_value / faceList.size();
	//		itor.second.postprocess = true;
	//	}
	//	else
	//	{
	//		itor.second.sdf_value = avg_sdf_value;
	//		itor.second.postprocess = true;
	//	}
	//}
	//BilateralSdfFilter(pMesh_);
	//
	//for (auto& itor : sdf_map)
	//{
	//	sdf_list3.push_back(itor.second.sdf_value);
	//}

	//std::sort(sdf_list3.begin(), sdf_list3.end());
	//float sdf_value_range = sdf_list3[sdf_list3.size() - 1] - sdf_list3[0];
	//for (auto& itor : sdf_map)
	//{
	//	itor.second.sdf_value = (itor.second.sdf_value - sdf_list3[0]) / sdf_value_range;
	//}
	//sdf_list3.clear();

	//for (auto& itor : sdf_map)
	//{
	//	sdf_list3.push_back(itor.second.sdf_value);
	//}

	///*node_->addDrawable(pMesh_);
	//pMesh_->setColorMode(CMPerVert);
	//for (auto& iter : sdf_map)
	//{
	//	if (iter.second.sdf_value <= 0.25)
	//	{
	//		iter.first->V(0)->C() = Color4b(255, 0, 0, 0);
	//		iter.first->V(1)->C() = Color4b(255, 0, 0, 0);
	//		iter.first->V(2)->C() = Color4b(255, 0, 0, 0);
	//	}
	//	else if (iter.second.sdf_value <= 0.5)
	//	{
	//		iter.first->V(0)->C() = Color4b(0, 255, 0, 0);
	//		iter.first->V(1)->C() = Color4b(0, 255, 0, 0);
	//		iter.first->V(2)->C() = Color4b(0, 255, 0, 0);
	//	}
	//	else if (iter.second.sdf_value <= 0.75)
	//	{
	//		iter.first->V(0)->C() = Color4b(0, 0, 255, 0);
	//		iter.first->V(1)->C() = Color4b(0, 0, 255, 0);
	//		iter.first->V(2)->C() = Color4b(0, 0, 255, 0);
	//	}
	//	else
	//	{
	//		iter.first->V(0)->C() = Color4b(0, 255, 255, 0);
	//		iter.first->V(1)->C() = Color4b(0, 255, 255, 0);
	//		iter.first->V(2)->C() = Color4b(0, 255, 255, 0);
	//	}
	//}*/


	return EXIT_SUCCESS;
}

//双边滤波
void CropScanRod::PImpl::BilateralSdfFilter(SmartMesh mesh)
{
	std::map<TFace*, SDF_property_map> new_sdf_map;
	float W = round(sqrt(mesh->numFace() / 2000)) + 1;
	float Rs = W / 2;
	std::vector<TFace*> ringFaceList;
	for (auto& itor : sdf_map)
	{
		ringFaceList.clear();
		mesh->getRingFace(itor.first, ringFaceList);
		float Rr = 0.0f;
		for (auto& iter : ringFaceList)
		{
			Rr += pow(sdf_map[iter].sdf_value - sdf_map[itor.first].sdf_value, 2);
		}
		Rr = sqrt(Rr / ringFaceList.size());
		float Ws = 0.0f, Wr = 0.0f;
		float weight_value = 0.0, new_sdf_value = 0.0;
		for (auto& iter : ringFaceList)
		{
			float centerLength = ((itor.first->V(0)->P() + itor.first->V(1)->P() + itor.first->V(2)->P()) / 3 - (iter->V(0)->P() + iter->V(1)->P() + iter->V(2)->P()) / 3).Length();
			Ws = exp(-pow(centerLength, 2) / (2 * pow(Rs, 2)));
			auto it = sdf_map.find(iter);
			if (it != sdf_map.end())
			{
				float sdfDiffValue = it->second.sdf_value - itor.second.sdf_value;
				Wr = exp(-pow(sdfDiffValue, 2) / (2 * pow(Rr, 2)));
			}
			new_sdf_value += Ws * Wr * itor.second.sdf_value;
			weight_value += Ws * Wr;
		}
		new_sdf_map[itor.first].sdf_value = new_sdf_value / weight_value;
	}
	sdf_map.clear();
	sdf_map.insert(new_sdf_map.begin(), new_sdf_map.end());
}

void CropScanRod::PImpl::sdf_values(const double boxRadius)
{

//	Vec3d face_centroid;
//	SDF_property_map face_property;
//	std::map<TFace*, std::vector<float>> ray_length;
//	MeshQuery::IntersectionList out_intersection_list;
//	sdf_map.clear();
//	mq_ = pMesh_->createTree(ITree::TT_QUERY)->asQuery();
//	for (auto& iter : pMesh_->face)
//	{
//		face_centroid = Vec3d::Construct((iter.V(0)->P() + iter.V(1)->P() + iter.V(2)->P()) / 3);
//		ray_.SetOrigin(face_centroid);
//		ray_.SetDirection(Vec3d::Construct(-iter.N()).normalized());
//		if (mq_->intersects(ray_, out_intersection_list))
//		{
//			/*Vec3d interP1 = out_intersection_list[0].intersection_point;
//			ray_.SetDirection(Vec3d::Construct(-iter.N()).normalized());
//			if (mq_->intersects(ray_, out_intersection_list))
//			{
//				Vec3d interP2 = out_intersection_list[0].intersection_point;
//				face_property.sdf_value = (interP1 - interP2).Length() - boxRadius * 2;
//				face_property.postprocess = true;
//				if(face_property.sdf_value<0)
//					face_property.postprocess = false;
//			}
//			else
//			{
//				face_property.postprocess = false;
//			}*/
//			Vec3d interP = out_intersection_list[0].intersection_point;
//			Vec3 interPlaneNor = Vec3::Construct(pMesh_->getFace(out_intersection_list[0].triangle_index)->N()).normalized();
//			float theta = asin((-iter.N().normalized() ^ interPlaneNor).Length());
//			face_property.sdf_value = (interP - face_centroid).Length() -2* boxRadius /*(1 / cos(theta) - 1) * boxRadius*/;
//			face_property.postprocess = true;
//			if (face_property.sdf_value < 0)
//				face_property.postprocess = false;
//			face_property.postprocess = true;
//		}
//		else
//		{
//			face_property.postprocess = false;
//		}
//		sdf_map[&iter] = face_property;
//	}
//
//	FaceList ringFaceList;
//	std::set<float> ringSdfList;
//	bool isMinSdf = false;
//	float avg_sdf_value = 0.0f,new_avg_sdf_value = 0.0f;
//	for (auto& iter : sdf_map)
//	{
//		if(iter.second.postprocess)
//			continue;
//		pMesh_->getRingFace(iter.first, ringFaceList);
//		for (auto& it : ringFaceList)
//		{
//			if (!sdf_map[it].postprocess)
//			{
//				isMinSdf = true;
//				continue;
//			}
//			ringSdfList.insert(sdf_map[it].sdf_value);
//			avg_sdf_value += sdf_map[it].sdf_value;
//		}
//		avg_sdf_value /= ringSdfList.size();
//		if(isMinSdf)
//		{
//			for (auto& it : ringFaceList)
//			{
//				if (!sdf_map[it].postprocess)
//				{
//					new_avg_sdf_value += avg_sdf_value;
//					sdf_map[it].sdf_value = avg_sdf_value;
//					sdf_map[it].postprocess = true;
//				}
//				new_avg_sdf_value += avg_sdf_value;
//			}
//			iter.second.sdf_value = new_avg_sdf_value / ringFaceList.size();
//			iter.second.postprocess = true;
//		}
//		else
//		{
//			iter.second.sdf_value = avg_sdf_value;
//			iter.second.postprocess = true;
//		}
//	}
//	for (auto& iter : sdf_map)
//	{
//		if (!iter.second.postprocess)
//			continue;
//	}
//
//	//BilateralSdfFilter(pMesh_);
//	node_->addDrawable(pMesh_);
//	pMesh_->setColorMode(CMPerVert);
//	std::vector<double> all_sdf;
//	for (auto& iter : sdf_map)
//	{
//		all_sdf.push_back(iter.second.sdf_value);
//		if (all_sdf.size() < floor(sdf_map.size() / 4))
//		{
//			iter.first->V(0)->C() = Color4b(255, 0, 0, 0);
//			iter.first->V(1)->C() = Color4b(255, 0, 0, 0);
//			iter.first->V(2)->C() = Color4b(255, 0, 0, 0);
//		}
//		else if (all_sdf.size() < floor(sdf_map.size() / 2))
//		{
//			iter.first->V(0)->C() = Color4b(0, 255, 0, 0);
//			iter.first->V(1)->C() = Color4b(0, 255, 0, 0);
//			iter.first->V(2)->C() = Color4b(0, 255, 0, 0);
//		}
//		else if (all_sdf.size() < floor(sdf_map.size() * 3 / 4))
//		{
//			iter.first->V(0)->C() = Color4b(0, 0, 255, 0);
//			iter.first->V(1)->C() = Color4b(0, 0, 255, 0);
//			iter.first->V(2)->C() = Color4b(0, 0, 255, 0);
//		}
//		else
//		{
//			iter.first->V(0)->C() = Color4b(0, 255, 255, 0);
//			iter.first->V(1)->C() = Color4b(0, 255, 255, 0);
//			iter.first->V(2)->C() = Color4b(0, 255, 255, 0);
//		}
//	}

#if 0
	sdf_map.clear();
	SDF_property_map cur_sdf;
	Vec3d face_centroid;
	std::vector<double> ray_length;
	mq_ = pMesh_->createTree(ITree::TT_QUERY)->asQuery();
	for (auto& iter : pMesh_->face)
	{
		face_centroid = Vec3d::Construct((iter.V(0)->P() + iter.V(1)->P() + iter.V(2)->P()) / 3);
		ray_.SetOrigin(face_centroid);
		//ray_.SetOrigin(face_centroid);
		//ray_.SetDirection(Vec3d::Construct(-iter.N()).normalized());
		//if (mq_->intersects(ray_, inter_))
		//{
		//	ray_length.push_back(inter_.dist-2* boxRadius);
		//}
		MeshQuery::IntersectionList inter_;
		float theta = 2 * M_PI / numbers_of_rays_;
		Vec3 zDir = -iter.N(), xDir = { 1.0f,0,0 }, yDir = { 0,1.0f,0 };
		XOYOZ(zDir, xDir, yDir);
		for (int i = 0; i < numbers_of_rays_; i++)
		{
			float angle = i * theta;
			ray_.SetDirection(Vec3d::Construct((xDir * sqrt(3) * std::rand() * sin(angle) + yDir * sqrt(3) * std::rand() * cos(angle)).normalized()));
			if (mq_->intersects(ray_, inter_))
			{
				ray_length.push_back(inter_[0].dist);
			}
		}
		//去除离群点
		/*sort(ray_length.begin(), ray_length.end());
		double median_value;
		if (ray_length.size() % 2 == 0)
		{
			median_value = (ray_length[ray_length.size() / 2 - 1] + ray_length[ray_length.size() / 2]) / 2;
		}
		else
		{
			median_value = ray_length[(ray_length.size() -1)/ 2];
		}*/
		if (ray_length.size() == 0)
		{
			cur_sdf.postprocess = false;
		}
		else
		{
			cur_sdf.postprocess = true;
			double avg_sdk = 0.0;
			for (auto iter : ray_length)
			{
				avg_sdk += iter;
			}
			cur_sdf.sdf_value = avg_sdk / ray_length.size();
		}
		sdf_map[&iter] = cur_sdf;
	}
	for (auto& itor : sdf_map)
	{
		if (itor.second.postprocess)
			continue;
		std::vector<TFace*> faceList;
		std::vector<double> ringSdfList;
		//遍历邻接面片
		pMesh_->getRingFace(itor.first, faceList);
		bool isMinValue = false;
		for (auto iter : faceList)
		{
			if (!sdf_map[iter].postprocess)
			{
				isMinValue = true;
				continue;
			}
			ringSdfList.push_back(sdf_map[iter].sdf_value);
		}
		if (isMinValue)
		{
			sort(ringSdfList.begin(), ringSdfList.end());
			itor.second.sdf_value = *ringSdfList.begin();
			itor.second.postprocess = true;
		}
		else
		{
			double avg_sdf = 0;
			for (auto ver : ringSdfList)
			{
				avg_sdf += ver;
			}
			avg_sdf = avg_sdf / ringSdfList.size();
			itor.second.sdf_value = avg_sdf;
			itor.second.postprocess = true;
		}
	}

	BilateralSdfFilter(pMesh_);
	node_->addDrawable(pMesh_);
	pMesh_->setColorMode(CMPerVert);
	std::vector<double> all_sdf;
	for (auto& iter : sdf_map)
	{
		all_sdf.push_back(iter.second.sdf_value);
		if (all_sdf.size() < floor(sdf_map.size() / 4))
		{
			iter.first->V(0)->C() = Color4b(255, 0, 0, 0);
			iter.first->V(1)->C() = Color4b(255, 0, 0, 0);
			iter.first->V(2)->C() = Color4b(255, 0, 0, 0);
		}
		else if (all_sdf.size() < floor(sdf_map.size() / 2))
		{
			iter.first->V(0)->C() = Color4b(0, 255, 0, 0);
			iter.first->V(1)->C() = Color4b(0, 255, 0, 0);
			iter.first->V(2)->C() = Color4b(0, 255, 0, 0);
		}
		else if (all_sdf.size() < floor(sdf_map.size()*3 / 4))
		{
			iter.first->V(0)->C() = Color4b(0, 0, 255, 0);
			iter.first->V(1)->C() = Color4b(0, 0, 255, 0);
			iter.first->V(2)->C() = Color4b(0, 0, 255, 0);
		}
		else 
		{
			iter.first->V(0)->C() = Color4b(0, 255, 255, 0);
			iter.first->V(1)->C() = Color4b(0, 255, 255, 0);
			iter.first->V(2)->C() = Color4b(0, 255, 255, 0);
		}
	}
#endif
	/*sort(all_sdf.begin(), all_sdf.end());
	for (auto& iter : sdf_map)
	{
		iter.second.sdf_value = (iter.second.sdf_value - all_sdf[0]) / (all_sdf[all_sdf.size() - 1] - all_sdf[0]);
	}

	node_->addDrawable(mesh);
	mesh->setColorMode(CMPerVert);*/
	/*for (auto& iter : mesh->face)
	{
		if (sdf_map[&iter].sdf_value < 0.25)
		{
			iter.V(0)->C() = Color4b(255, 0, 0, 0);
			iter.V(1)->C() = Color4b(255, 0, 0, 0);
			iter.V(2)->C() = Color4b(255, 0, 0, 0);
		}
		else if (sdf_map[&iter].sdf_value < 0.5)
		{
			iter.V(0)->C() = Color4b(0, 255, 0, 0);
			iter.V(1)->C() = Color4b(0, 255, 0, 0);
			iter.V(2)->C() = Color4b(0, 255, 0, 0);
		}
		else if(sdf_map[&iter].sdf_value < 0.75)
		{
			iter.V(0)->C() = Color4b(0, 0, 255, 0);
			iter.V(1)->C() = Color4b(0, 0, 255, 0);
			iter.V(2)->C() = Color4b(0, 0, 255, 0);
		}
		else
		{
			iter.V(0)->C() = Color4b(0, 255, 255, 0);
			iter.V(1)->C() = Color4b(0, 255, 255, 0);
			iter.V(2)->C() = Color4b(0, 255, 255, 0);
		}
	}*/
}
#if 0
//双边滤波
void CropScanRod::PImpl::BilateralSdfFilter(SmartMesh mesh)
{
	std::map<TFace*, SDF_property_map> new_sdf_map;
	float W = round(sqrt(mesh->numFace() / 2000)) + 1;
	float Rs = W / 2;
	std::vector<TFace*> ringFaceList;
	for (auto& itor : sdf_map)
	{
		ringFaceList.clear();
		mesh->getRingFace(itor.first, ringFaceList);
		float Rr = 0.0f;
		for (auto& iter : ringFaceList)
		{
			Rr += pow(sdf_map[iter].sdf_value - sdf_map[itor.first].sdf_value, 2);
		}
		Rr = sqrt(Rr/ringFaceList.size());
		float Ws = 0.0f, Wr = 0.0f;
		float weight_value = 0.0, new_sdf_value = 0.0;
		for (auto& iter : ringFaceList)
		{
			float centerLength = ((itor.first->V(0)->P() + itor.first->V(1)->P() + itor.first->V(2)->P()) / 3 - (iter->V(0)->P() + iter->V(1)->P() + iter->V(2)->P()) / 3).Length();
			Ws = exp(-pow(centerLength, 2) / (2 * pow(Rs, 2)));
			auto it = sdf_map.find(iter);
			if (it != sdf_map.end())
			{
				float sdfDiffValue = it->second.sdf_value - itor.second.sdf_value;
				Wr = exp(-pow(sdfDiffValue, 2) / (2 * pow(Rr, 2)));
			}
			new_sdf_value += Ws * Wr * itor.second.sdf_value;
			weight_value += Ws * Wr;
		}
		new_sdf_map[itor.first].sdf_value = new_sdf_value / weight_value;
	}
	sdf_map.clear();
	sdf_map.insert(new_sdf_map.begin(), new_sdf_map.end());
}

void CropScanRod::PImpl::SoftCluster()
{
	//选取k个种子点
	shared_ptr<ObbFilter> obbObject = std::make_shared<ObbFilter>(pMesh_.get());
	float a = obbObject->getWidth();
	float b = obbObject->getHeight();
	float c = obbObject->getDepth();
	Vec3 aDir = obbObject->getAxes().main;
	Vec3 bDir = obbObject->getAxes().second;
	Vec3 cDir = obbObject->getAxes().third;
	std::pair<float, Vec3f> maxAxis = a > b ? (a > c ? std::make_pair(a, aDir) : std::make_pair(c, cDir)) : (b > c ? std::make_pair(b, bDir) : std::make_pair(c, cDir));
	Vec3 perDir = maxAxis.second * maxAxis.first / typeNum_;

	std::map<int,TFace*> random_seeds;
	Vec3 curPts = obbObject->getCenter() - maxAxis.second * maxAxis.first / 2 + perDir/2;
	//auto copy_mesh = Mesh::copyMesh(*pMesh_);
	ann_ = pMesh_->createTree(ITree::TT_ANN)->asAnn();
	std::vector<int> out_index;
	std::vector<float> out_dist;
	for (int i = 0; i < typeNum_; i++)
	{
		out_index.clear();
		out_dist.clear();
		Vec3 queryPt = curPts + perDir * i;
		ann_->query(queryPt, 1, out_index, out_dist);
		TVertex* vertIter = pMesh_->getVertex(out_index[0]);
		VFIterator fv(vertIter);
		random_seeds[i] = fv.F();
		node_->addDrawable(SmartDrawable(new PointGeometry(vertIter->P(), Color4f::Red, 5.0f)));
	}

	if (random_seeds.size() < typeNum_)
		return;

	float minDist = FLT_MAX;
	float before_SSE = FLT_MAX;
	float before_SC = FLT_MIN;
	float diff_value = 0.0f;
	bool bCluster_ = true;
	while (bCluster_)
	{
		for (auto& iter : sdf_map)
		{
			for (auto& itor : random_seeds)
			{
				diff_value = abs(iter.second.sdf_value - sdf_map[itor.second].sdf_value);
				if (diff_value < minDist)
				{
					minDist = diff_value;
					iter.second.output_cluster_ids = true;
					iter.second.number_of_clusters = itor.first;
				}
			}
		}


		std::map<int, std::vector<TFace*>> clusterFaceList;
		std::map<int, float> newClusterSdfList;//每个聚类对应的平均sdf值
		for (auto& iter : sdf_map)
		{
			clusterFaceList[iter.second.number_of_clusters].push_back(iter.first);
		}
		for (auto& iter : clusterFaceList)
		{
			float new_cluster_sdf = 0.0f;
			for (auto& itor : iter.second)
			{
				new_cluster_sdf += sdf_map[itor].sdf_value;
			}
			new_cluster_sdf /= iter.second.size();

			newClusterSdfList[iter.first] = new_cluster_sdf;
		}

		float SSE = 0.0f, SC = 0.0f; //SSE、轮廓系数
		float cur_cluster_avg_dist = 0.0f;
		float near_cluster_avg_dist = 0.0f;

		for (auto& iter : clusterFaceList)
		{
			for (auto& itor : iter.second)
			{
				//样本与同一簇的其他样本点的平均距离
				for (auto& it : iter.second)
				{
					if (itor == it)
						continue;
					cur_cluster_avg_dist += abs(sdf_map[itor].sdf_value - sdf_map[it].sdf_value);
					SSE += pow((sdf_map[itor].sdf_value - sdf_map[it].sdf_value),2);//每个样本点到所在簇的中心的误差平方和
				}
				minDist = FLT_MAX;
				int near_cluster = 0;
				for (auto& rs : random_seeds)
				{
					if (sdf_map[itor].number_of_clusters != rs.first)
					{
						diff_value = abs(sdf_map[itor].sdf_value - sdf_map[rs.second].sdf_value);
						if (diff_value < minDist)
						{
							minDist = diff_value;
							near_cluster = rs.first;
						}
					}
				}
				cur_cluster_avg_dist /= iter.second.size();
				//样本与距离最近簇中所有样本点的平均距离
				for (auto& cf : clusterFaceList[near_cluster])
				{
					if (sdf_map[itor].number_of_clusters == near_cluster)
						continue;
					near_cluster_avg_dist += abs(sdf_map[itor].sdf_value - sdf_map[cf].sdf_value);
				}
				near_cluster_avg_dist /= clusterFaceList[near_cluster].size();

				//当前样本点的轮廓系数
				SC += (near_cluster_avg_dist - cur_cluster_avg_dist) / std::max(cur_cluster_avg_dist, near_cluster_avg_dist);
			}

		}
		SC /= sdf_map.size();
		
		//计算新的中心
		minDist = FLT_MAX;
		random_seeds.clear();
		for (auto& itor : newClusterSdfList)
		{
			for (auto& iter : sdf_map)
			{
				diff_value = abs(iter.second.sdf_value - itor.second);
				SSE += pow(diff_value, 2);
				if (diff_value < minDist)
				{
					minDist = diff_value;
					random_seeds[itor.first] = iter.first;
				}
			}
		}
		
		if (SC > before_SC || SSE < before_SSE)
		{
			before_SC = SC;
			before_SSE = SSE;
		}
		else
		{
			bCluster_ = false;
		}
	}

	/*std::map<TFace*, float> probability_value_list;
	for (auto& iter : sdf_map)
	{
		probability_value_list[iter.first] = DistributionProbability(iter.second.number_of_clusters);
	}*/
}
//概率分布
float CropScanRod::PImpl::DistributionProbability(std::size_t cluster)
{
	//期望最大 Expectation Maximization Algorithm
	//最大似然估计可得正态分布的模型参数在MLE下的估计值
	//mu = (x1+x2+....+xn)/n;
	//sigma = [(x1-mu)^2+(x2-mu)^2+....+(xn-mu)^2]/n
	float mu = 0.0f, sigma = 0.0f;
	for (auto& iter : sdf_map)
	{
		mu += iter.second.number_of_clusters;
	}
	mu /= sdf_map.size();
	for (auto& iter : sdf_map)
	{
		sigma += pow((cluster - mu), 2);
	}
	sigma = sqrt(sigma / sdf_map.size());
	//概率分布(高斯分布)函数： P(x) = 1/sqrt(2*pi*sigma^2) * e^(-(x-mu)^2/(2*sigma^2))
	float a = 1 / sqrt(2 * M_PI * pow(sigma, 2));
	float b = exp(-pow(cluster - mu, 2) / (2 * pow(sigma, 2)));
	return a * b;
}

double CropScanRod::PImpl::alpha_expansion_graphcut(float threshold, double smoothing_lambda)
{
	Graph graph;
	float maxflow = FLT_MAX, curflow;
	std::map<int, std::vector<TFace*>> clusterFaceList;

	std::map<int, std::vector<TFace*>> result_label;
	std::set<int> alpha_list;
	for (auto& iter:sdf_map)
	{
		clusterFaceList[iter.second.number_of_clusters].push_back(iter.first);
		alpha_list.insert(iter.second.number_of_clusters);
	}
	float forward_weight_tlink, backward_weight_tlink, forward_weight_nlink, backward_weight_nlink;
	for (auto alpha : alpha_list)
	{
		SelectFilter::clearSelect(pMesh_.get());
		for (auto& iter : clusterFaceList)
		{
			float probability_value_p = DistributionProbability(iter.first);//概率分布
			float probability_value_alpha = DistributionProbability(alpha);
			forward_weight_tlink = -log(std::max(probability_value_alpha, threshold));
			if (iter.first == alpha)
			{
				backward_weight_tlink = FLT_MAX;
			}
			else
			{
				backward_weight_tlink = -log(std::max(probability_value_p, threshold));
			}
			// 添加t_link权重
			for (auto itor : iter.second)
			{
				if(itor->IsS())
					continue;
				graph.addTLink(itor->id(), forward_weight_tlink, backward_weight_tlink);
				itor->SetS();
				// 添加n_link权重
				FaceList ringFaceList;
				pMesh_->getRingFace(itor, ringFaceList);
				for (auto it : ringFaceList)
				{
					if (sdf_map[it].number_of_clusters == iter.first)
					{
						forward_weight_nlink = 0;
					}
					else
					{
						float cone_angle = sqrt(1 - pow(itor->N() * it->N(), 2));//二面角的正弦
						forward_weight_nlink = cone_angle >= 0 ? -log(asin(cone_angle) / M_PI) : -log(0.1 * asin(cone_angle) / M_PI);
						int auxiliary_node_id = graph.addVertex();
						graph.addTLink(auxiliary_node_id, 0, forward_weight_nlink);
						if (iter.first == alpha)
							graph.addNLink(itor->id(), auxiliary_node_id, 0, 0);
						else
							graph.addNLink(itor->id(), auxiliary_node_id, 0, forward_weight_nlink);

						if (sdf_map[it].number_of_clusters == alpha)
							graph.addNLink(auxiliary_node_id, it->id(), 0, 0);
						else
							graph.addNLink(auxiliary_node_id, it->id(),0, forward_weight_nlink);
					}
				}
			}
			curflow = graph.maxFlow();
			for (auto itor : iter.second)
			{
				if (graph.inSourceSegment(itor->id())==1)
				{
					result_label[alpha].push_back(itor);
				}
			}
		}
	}
	return curflow;
}


void CropScanRod::PImpl::EnergyFunction(float threshold, double smoothing_lambda)
{
	std::map<int, std::vector<TFace*>> clusterFaceList;
	for (auto& iter : sdf_map)
	{
		clusterFaceList[iter.second.number_of_clusters].push_back(iter.first);
	}
	float before_expansion_energy = 0.0f, after_expansion_energy = 0.0f;
	bool success = false;
	std::map<TFace*, SDF_property_map> after_sdf_map;
	copy(sdf_map.begin(), sdf_map.end(), inserter(after_sdf_map, after_sdf_map.begin()));
	while (!success)
	{
		float e_data, e_smooth;//e1标签函数，e2平滑程度
		float cone_angle = 0.0f;
		std::vector<TFace*> ringFace;
		for (auto& iter : clusterFaceList)
		{
			float probability_value = DistributionProbability(iter.first);//概率分布
			for (auto& itor : iter.second)
			{
				e_data = -log(std::max(probability_value, threshold));
				before_expansion_energy += e_data;
				after_expansion_energy += e_data;
				pMesh_->getRingFace(iter.first, ringFace);
				for (auto& it : ringFace)
				{
					if (sdf_map[itor].number_of_clusters == sdf_map[it].number_of_clusters)
					{
						e_smooth = 0;
					}
					else
					{
						cone_angle = sqrt(1 - pow(itor->N() * it->N(), 2));//二面角的余弦
						e_smooth = cone_angle >= 0 ? -log(cone_angle / M_PI) : -log(0.1 * cone_angle / M_PI);
					}
					after_expansion_energy += smoothing_lambda * e_smooth;
				}
			}
		}
		if (after_expansion_energy < before_expansion_energy)
		{
			success = true;
			sdf_map.clear();
			copy(after_sdf_map.begin(), after_sdf_map.end(), inserter(sdf_map, sdf_map.begin()));
		}
	}
	//The energy function minimized using alpha-expansion graph cut algorithm （图割之 Alpha-expansion）
}
#endif

void CropScanRod::PImpl::XOYOZ(Vec3& zDir, Vec3& xDir, Vec3& yDir)
{
	Matrix rot;
	float angle = Angle(zDir, { 0.0f,0.0f,1.0f });
	Vec3 rotDir = Vec3(0.0f, 0.0f, 1.0f) ^ zDir;
	rot.SetRotateRad(angle, rotDir);
	//向量a,b
	xDir = rot * xDir;
	yDir = rot * yDir;
	xDir = xDir.normalized();
	yDir = yDir.normalized();
}


