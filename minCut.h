#ifndef SURFACE_MESH_SEGMENTATION_H
#define SURFACE_MESH_SEGMENTATION_H
#include<Mesh.h>
using namespace std;
typedef unsigned char uchar;
/////////////////////////////////////////////////////////
//This section shows how to use the library to compute
//a minimum cut on the following graph :
//                         SOURCE
//                        /       \
//                      1/         \2
//                      /      3    \
//                    node0 -----> node1
//                      |   <-----   |
//                      |      4     |
//                      \            /
//                      5\          /6
//                        \        /
//                           SINK
///////////////////////////////////////////////////
class Graph
{
public:
	enum termtype {
		SOURCE = 0,
		SINK = 1,
	};


	Graph();
	Graph(int VertexCount, int EdgeCount);
	~Graph();
	void create(int VertexCount, int EdgeCount);
	int addVertex();
	void addNLink(int i, int j, float w1, float w2);
	void addTLink(int i, float sourceW, float sinkW);
	float maxFlow();
	termtype inSourceSegment(int i);

private:
	class pointNode
	{
	public:
		pointNode* next;
		float weight;
		int time;
		int firstEdge;//首个相邻边
		int dist;//路径长度
		int parentNode;// describes node's parent
		uchar flowDir;//0为源节点，1为汇节点
	};
	class edgeNode
	{
	public:
		int nextNodeID;//边指向的结点
		int nextEdgeID;//改边的顶点的下一条边
		float weight;
	};

	std::vector<pointNode> pointNodeList;
	std::vector<edgeNode> edgeNodeList;
	float flow;//图的流量
};
Graph::Graph()
{
	flow = 0;
}

Graph::Graph(int VertexCount, int EdgeCount)
{
	create(VertexCount, EdgeCount);
}

Graph::~Graph()
{
}

//给图的节点容器和边容器分配内存空间
void Graph::create(int VertexCount, int EdgeCount)
{
	pointNodeList.reserve(VertexCount);
	edgeNodeList.reserve(EdgeCount + 2); 
	flow = 0;
}
//添加空节点
int Graph::addVertex()
{ 
	pointNode pn;
	memset(&pn, 0, sizeof(pointNode));
	pointNodeList.push_back(pn);
	return pointNodeList.size() - 1;
}
/*
添加点之间的n - link;
参数说明：
int-- - i: 弧头结点编号
int-- - j : 弧尾结点编号
float-- - w1 : 正向弧权值
float-- - w2 : 逆向弧权值
返回值：无
*/
void Graph::addNLink(int i, int j, float w1, float w2)
{
	assert(i >= 0 && i < pointNodeList.size());
	assert(j >= 0 && j < pointNodeList.size());
	assert(w1 >= 0 && w2 >= 0);
	assert(i != j);
	edgeNode fromEdge, toEdge;//正反向弧
	fromEdge.nextNodeID = j;//正向弧指向的节点
	fromEdge.nextEdgeID = pointNodeList[i].firstEdge;//
	fromEdge.weight = w1;
	pointNodeList[i].firstEdge = edgeNodeList.size();//修改节点i的第一个弧当做正向弧
	edgeNodeList.push_back(fromEdge);

	toEdge.nextNodeID = i;
	toEdge.nextEdgeID = pointNodeList[j].firstEdge;
	toEdge.weight = w2;
	pointNodeList[j].firstEdge = edgeNodeList.size();
	edgeNodeList.push_back(toEdge);
}


/*
添加节点到顶点的边的t-link
参数说明：
int-- - i: 结点编号
float-- - sourceW : 正向弧权值
float-- - sinkW : 逆向弧权值
返回值：无
*/
void Graph::addTLink(int i, float sourceW, float sinkW)
{
	assert(i >= 0 && i < pointNodeList.size());
	float weight = pointNodeList[i].weight;
	if (weight > 0)
		sourceW += weight;
	else
		sinkW -= weight;
	flow += (sourceW < sinkW) ? sourceW : sinkW;
	pointNodeList[i].weight = sourceW - sinkW;
}

//最大流
float Graph::maxFlow()
{
	const int TERMINAL = -1;//终点
	const int ORPHAN = -2;//孤点
	pointNode stub, *nilNode = &stub, *first = nilNode, *last = nilNode;//stub哨兵结点
	int curTime = 0;
	stub.next = nilNode;//首节点指向自己
	pointNode* pointPtr = &pointNodeList[0];//结点指针
	edgeNode* edgePtr = &edgeNodeList[0];//弧指针
	std::vector<pointNode*> orphanNodeList;//孤立点合集
	//遍历结点，初始化活动结点
	for (int i = 0; i < pointNodeList.size(); i++)
	{
		pointNode* pn = pointPtr + i;
		pn->time = 0;
		if (pn->weight != 0)
		{
			last = last->next = pn;
			pn->dist = 1;
			pn->parentNode = TERMINAL; // 标注其双亲为终端结点
			pn->flowDir = pn->weight < 0;
		}
		else
		{
			pn->parentNode = 0;//孤结点
		}
	}
	first = first->next;
	last->next = nilNode;//stub哨兵结点放置队尾
	nilNode->next = 0;

	
	while (true)
	{
		pointNode *curNode, *adjcencyNode;
		float minWeight, curWeight;
		uchar curflowDir;//正向0,反向1
		int e0 = -1, ei = 0, ej = 0;
		//----------------------------------source->sink第一条路径-----------------------------------------------
		while (first != nilNode)
		{
			curNode = first;
			if (curNode->parentNode)
			{
				curflowDir = curNode->flowDir;
				for (ei = curNode->firstEdge; ei != 0; ei = edgePtr[ei].nextEdgeID)
				{
					if (edgeNodeList[ei ^ curflowDir].weight == 0)
						continue;
					adjcencyNode = pointPtr + edgePtr[ei].nextNodeID;
					if (!adjcencyNode->parentNode)
					{
						adjcencyNode->flowDir = curflowDir;
						adjcencyNode->parentNode = ei ^ 1;
						adjcencyNode->time = curNode->time;
						adjcencyNode->dist += 1;
						if (!adjcencyNode->next)
						{
							adjcencyNode->next = nilNode;
							last = last->next = adjcencyNode;
						}
						continue;
					}

					if (adjcencyNode->flowDir != curflowDir)
					{
						e0 = ei ^ curflowDir;
						break;
					}

					if (adjcencyNode->dist > curNode->dist && adjcencyNode->time <= curNode->time)
					{
						adjcencyNode->parentNode = ei ^ 1;
						adjcencyNode->time = curNode->time;
						adjcencyNode->dist += 1;
					}
				}
				if (e0 > 0)
				{
					break;
				}

			}
			first = first->next;
			curNode->next = 0;
		}
		if (e0 <= 0)
		{
			break;
		}


		//----------------------------------流量统计与拆分-----------------------------------------------
		//第一节： 查找路径中的最小权值
		minWeight = edgePtr[e0].weight;
		assert(minWeight > 0);
		//遍历整条路径，从当前节点开始，向前回溯s树 向后回溯t树
		// 2次遍历， k=1: 回溯s树， k=0: 回溯t树
		for (int k = 1; k >= 0; k--)
		{
			for (curNode = pointPtr + edgePtr[e0 ^ k].nextNodeID;; curNode = pointPtr + edgePtr[ei].nextNodeID)
			{
				if ((ei = curNode->parentNode) < 0)
					break;
				curWeight = edgePtr[ei ^ k].weight;
				minWeight = min(minWeight, curWeight);
				assert(minWeight > 0);
			}
			curWeight = fabs(curNode->weight);
			minWeight = min(minWeight, curWeight);
			assert(minWeight > 0);
		}

		/*第二节：修改当前路径中的所有的weight权值
		任何时候s和t树的结点都只有一条边使其连接到树中，当这条弧权值减少为0则此结点从树中断开，
		若其无子结点，则成为孤立点，若其拥有子结点，则独立为森林，但是ei的子结点还不知道他们被孤立了！
		*/
		edgePtr[e0].weight -= minWeight;
		edgePtr[e0 ^ 1].weight += minWeight; 
		flow += minWeight;


		// k = 1: source tree, k = 0: destination tree
		for (int k = 1; k >= 0; k--)
		{
			for (curNode = pointPtr + edgePtr[e0 ^ k].nextNodeID;; curNode = pointPtr + edgePtr[ei].nextNodeID)
			{
				if ((ei = curNode->parentNode) < 0)
					break;
				edgePtr[ei ^ (k ^ 1)].weight += minWeight;
				if ((edgePtr[ei ^ k].weight -= minWeight) == 0)
				{
					orphanNodeList.push_back(curNode);
					curNode->parentNode = ORPHAN;
				}
			}
			curNode->weight = curWeight + minWeight*(1 - k * 2);
			if (curNode->weight == 0)
			{
				orphanNodeList.push_back(curNode);
				curNode->parentNode = ORPHAN;
			}
		}

		//---------------------------- 第三阶段: 树的重构 寻找新的父节点，恢复搜索树 -----------------------------//
		//为每个孤点p找到新的合法父节点q
		//p和q必须为一颗搜索树，两节点间的流量大于0 节点q的根节点为源点s或汇点t
		curTime++;
		while (!orphanNodeList.empty())
		{
			pointNode* curNode = orphanNodeList.back();
			orphanNodeList.pop_back();

			int d, minDist = INT_MAX;
			e0 = 0;
			curflowDir = curNode->flowDir;

			//遍历当前结点的相邻点
			for (ei = curNode->firstEdge; ei != 0; ei = edgePtr[ei].nextEdgeID)
			{
				if (edgePtr[ei ^ (curflowDir ^ 1)].weight == 0)
					continue;
				adjcencyNode = pointPtr + edgePtr[ei].nextNodeID;
				if (adjcencyNode->flowDir != curflowDir||adjcencyNode->parentNode==0)
					continue;

				//计算当前路径长度
				for (d = 0;;)
				{
					if (adjcencyNode->time == curTime)
					{
						d += adjcencyNode->dist;
						break;
					}
					ej = adjcencyNode->parentNode;
					d++;
					if (ej < 0)
					{
						if (ej == ORPHAN)
							d = INT_MAX - 1;
						else
						{
							adjcencyNode->time = curTime;
							adjcencyNode->dist = 1;
						}
						break;
					}
					adjcencyNode = pointPtr + edgePtr[ej].nextNodeID;
				}

				//update the distance
				if (++d < INT_MAX)
				{
					if (d < minDist)
					{
						minDist = d;
						e0 = ei;
					}
					for (adjcencyNode = pointPtr + edgePtr[ei].nextNodeID; adjcencyNode->time != curTime; adjcencyNode = pointPtr + edgePtr[adjcencyNode->parentNode].nextNodeID)
					{
						adjcencyNode->time = curTime;
						adjcencyNode->dist = --d;
					}
				}
			}
				
			if ((curNode->parentNode = e0) > 0)
			{
				curNode->time = curTime;
				curNode->dist = minDist;
				continue;
			}


			//没有找到该节点的合法父节点，节点变为自由节点，其子节点变为孤点
			curNode->time = 0;
			for(ei = curNode->firstEdge; ei!= 0; ei = edgePtr[ei].nextEdgeID)
			{
				adjcencyNode = pointPtr + edgePtr[ei].nextNodeID;
				ej = adjcencyNode->parentNode;
				if(adjcencyNode->flowDir!=curflowDir||!ej)
					continue;
				if (edgePtr[ei ^ (curflowDir ^ 1)].weight && !adjcencyNode->next)
				{
					adjcencyNode->next = nilNode;
					last = last->next = adjcencyNode;
				}
				if (ej > 0 && pointPtr + edgePtr[ej].nextNodeID == curNode)
				{
					orphanNodeList.push_back(adjcencyNode);
					adjcencyNode->parentNode = ORPHAN;
				}
			}
		}
	}
	return flow;
}

Graph::termtype Graph::inSourceSegment(int i)
{
	assert(i >= 0 && i < pointNodeList.size());
	if (pointNodeList[i].flowDir == 0)
		return SINK;
	else
		return SOURCE;
}
#endif 
