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
		int firstEdge;//�׸����ڱ�
		int dist;//·������
		int parentNode;// describes node's parent
		uchar flowDir;//0ΪԴ�ڵ㣬1Ϊ��ڵ�
	};
	class edgeNode
	{
	public:
		int nextNodeID;//��ָ��Ľ��
		int nextEdgeID;//�ıߵĶ������һ����
		float weight;
	};

	std::vector<pointNode> pointNodeList;
	std::vector<edgeNode> edgeNodeList;
	float flow;//ͼ������
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

//��ͼ�Ľڵ������ͱ����������ڴ�ռ�
void Graph::create(int VertexCount, int EdgeCount)
{
	pointNodeList.reserve(VertexCount);
	edgeNodeList.reserve(EdgeCount + 2); 
	flow = 0;
}
//��ӿսڵ�
int Graph::addVertex()
{ 
	pointNode pn;
	memset(&pn, 0, sizeof(pointNode));
	pointNodeList.push_back(pn);
	return pointNodeList.size() - 1;
}
/*
��ӵ�֮���n - link;
����˵����
int-- - i: ��ͷ�����
int-- - j : ��β�����
float-- - w1 : ����Ȩֵ
float-- - w2 : ����Ȩֵ
����ֵ����
*/
void Graph::addNLink(int i, int j, float w1, float w2)
{
	assert(i >= 0 && i < pointNodeList.size());
	assert(j >= 0 && j < pointNodeList.size());
	assert(w1 >= 0 && w2 >= 0);
	assert(i != j);
	edgeNode fromEdge, toEdge;//������
	fromEdge.nextNodeID = j;//����ָ��Ľڵ�
	fromEdge.nextEdgeID = pointNodeList[i].firstEdge;//
	fromEdge.weight = w1;
	pointNodeList[i].firstEdge = edgeNodeList.size();//�޸Ľڵ�i�ĵ�һ������������
	edgeNodeList.push_back(fromEdge);

	toEdge.nextNodeID = i;
	toEdge.nextEdgeID = pointNodeList[j].firstEdge;
	toEdge.weight = w2;
	pointNodeList[j].firstEdge = edgeNodeList.size();
	edgeNodeList.push_back(toEdge);
}


/*
��ӽڵ㵽����ıߵ�t-link
����˵����
int-- - i: �����
float-- - sourceW : ����Ȩֵ
float-- - sinkW : ����Ȩֵ
����ֵ����
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

//�����
float Graph::maxFlow()
{
	const int TERMINAL = -1;//�յ�
	const int ORPHAN = -2;//�µ�
	pointNode stub, *nilNode = &stub, *first = nilNode, *last = nilNode;//stub�ڱ����
	int curTime = 0;
	stub.next = nilNode;//�׽ڵ�ָ���Լ�
	pointNode* pointPtr = &pointNodeList[0];//���ָ��
	edgeNode* edgePtr = &edgeNodeList[0];//��ָ��
	std::vector<pointNode*> orphanNodeList;//������ϼ�
	//������㣬��ʼ������
	for (int i = 0; i < pointNodeList.size(); i++)
	{
		pointNode* pn = pointPtr + i;
		pn->time = 0;
		if (pn->weight != 0)
		{
			last = last->next = pn;
			pn->dist = 1;
			pn->parentNode = TERMINAL; // ��ע��˫��Ϊ�ն˽��
			pn->flowDir = pn->weight < 0;
		}
		else
		{
			pn->parentNode = 0;//�½��
		}
	}
	first = first->next;
	last->next = nilNode;//stub�ڱ������ö�β
	nilNode->next = 0;

	
	while (true)
	{
		pointNode *curNode, *adjcencyNode;
		float minWeight, curWeight;
		uchar curflowDir;//����0,����1
		int e0 = -1, ei = 0, ej = 0;
		//----------------------------------source->sink��һ��·��-----------------------------------------------
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


		//----------------------------------����ͳ������-----------------------------------------------
		//��һ�ڣ� ����·���е���СȨֵ
		minWeight = edgePtr[e0].weight;
		assert(minWeight > 0);
		//��������·�����ӵ�ǰ�ڵ㿪ʼ����ǰ����s�� ������t��
		// 2�α����� k=1: ����s���� k=0: ����t��
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

		/*�ڶ��ڣ��޸ĵ�ǰ·���е����е�weightȨֵ
		�κ�ʱ��s��t���Ľ�㶼ֻ��һ����ʹ�����ӵ����У���������Ȩֵ����Ϊ0��˽������жϿ���
		�������ӽ�㣬���Ϊ�����㣬����ӵ���ӽ�㣬�����Ϊɭ�֣�����ei���ӽ�㻹��֪�����Ǳ������ˣ�
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

		//---------------------------- �����׶�: �����ع� Ѱ���µĸ��ڵ㣬�ָ������� -----------------------------//
		//Ϊÿ���µ�p�ҵ��µĺϷ����ڵ�q
		//p��q����Ϊһ�������������ڵ�����������0 �ڵ�q�ĸ��ڵ�ΪԴ��s����t
		curTime++;
		while (!orphanNodeList.empty())
		{
			pointNode* curNode = orphanNodeList.back();
			orphanNodeList.pop_back();

			int d, minDist = INT_MAX;
			e0 = 0;
			curflowDir = curNode->flowDir;

			//������ǰ�������ڵ�
			for (ei = curNode->firstEdge; ei != 0; ei = edgePtr[ei].nextEdgeID)
			{
				if (edgePtr[ei ^ (curflowDir ^ 1)].weight == 0)
					continue;
				adjcencyNode = pointPtr + edgePtr[ei].nextNodeID;
				if (adjcencyNode->flowDir != curflowDir||adjcencyNode->parentNode==0)
					continue;

				//���㵱ǰ·������
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


			//û���ҵ��ýڵ�ĺϷ����ڵ㣬�ڵ��Ϊ���ɽڵ㣬���ӽڵ��Ϊ�µ�
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
