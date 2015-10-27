#ifndef GRAPHWINDOW_H
#define GRAPHWINDOW_H
 
#include <QVTKWidget.h>
#include <QtGui/QMainWindow>
#include <vtkTable.h>
#include <vtkTableToGraph.h>
#include <vtkGraphLayoutView.h>
#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>
#include <vtkLookupTable.h>
#include <vtkObject.h>
#include <vtkPoints.h>
#include "ObjectSelection.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include "ftkGUI/ColorMap.h"
#include <vtkCornerAnnotation.h>

#ifndef MYPOINT
#define MYPOINT
typedef struct Point
{
	Point(){};
	Point( double ix, double iy)
	{
		x = ix;
		y = iy;
	};

	double x;
	double y;
}Point;
#endif

class GraphWindow : public QMainWindow
{
    Q_OBJECT;

public:
	GraphWindow(QWidget * parent = 0);
	~GraphWindow();
	void setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels = NULL, std::vector<int> *indexCluster = NULL, ObjectSelection * sels2 = NULL);
	vtkSmartPointer<vtkTable> getModels();
	void SetGraphTable(vtkSmartPointer<vtkTable> table);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::string xCol, std::string yCol, std::string zCol);
	void SetTreeTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::vector<double> *colorVec = NULL, 
					std::vector<double> *disVec = NULL, std::set<long int>* colSels = NULL, QString filename = "");
	void AdjustedLayout(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::vector<int> *treeOrder = NULL, std::vector<double> *colorVec = NULL, std::vector<double> *disVec = NULL);
	void ShowGraphWindow();
	ObjectSelection * GetSelection();
	void GetTrendTreeOrder(std::vector<long int> &order);
	void ColorTreeAccordingToFeatures(vnl_vector<double> &feature, const char *featureName);
	static void GetTreeNodeBetweenDistance(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, vnl_matrix<double> &disMat);
	void GetLongestPathListByLevel(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, int level, std::vector<std::list<int>> &longestPath);
	void GetLongestPathListCluster(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::vector<std::list<int>> &longestCluster);

protected:
	void SetSelectedIds(std::set<int>& IDs);
	void SetSelectedIds2();
	void UpdataLookupTable( std::set<long int>& IDs);
	void CalculateCoordinates(std::vector< std::list<int>> &adjList, std::vector<Point>& pointList);
	void find(vnl_vector<int>& vec, int val, std::vector<int>& equal, std::vector<int>& nonequal);
	void getBackBones(vnl_matrix< int>& shortest_hop, std::vector< int>& branchnodes, std::vector< int>& chains);
	void GetElementsIndexInMatrix( vnl_matrix<double>& nodePos, vnl_matrix<double>& repel_mat, vnl_vector< int>& tag);
	void SortChainList( std::vector<int>& backbones, std::vector< std::pair<int, std::vector<int> > >& chainList);
	int IsExist(std::vector<int>& vec, int value);
	double Median( vnl_vector<double> vec);
	void GetOrder(int node, std::vector<long int> &order);
	void UpdateCoordinatesByEdgeWeights(std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList);
	void UpdateChainPointList(int attachnode, std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList);
	Point GetNewPointFromOldPoint( Point &oldPointFirst, Point &oldPointSecond, Point &newPointFirst, double weight);
	double GetEdgeWeight(vnl_matrix<double>& vertexList, long firstIndex, long secondIndex);
	virtual void closeEvent(QCloseEvent *event);
	void SetUserDefineTrend(int nodeID);
	void SetTrendStartTag(bool bstart);
	void UpdateTrendPath();
	void GetTrendPath(vnl_matrix<int> &hopMat, int startNode, int endNode, std::vector< int> &path);
	void ResetLookupTable(vtkSmartPointer<vtkLookupTable> lookuptable, double* color);
	void RestoreLookupTable();

	bool FindCycle(vnl_matrix<unsigned char> &adjacent_matrix, std::vector< std::list<int> > &cycleList, std::vector<int> &seperatePts);
	void SearchCycles(vnl_matrix<unsigned char> &adjacent_matrix, std::list<int> &cycleNode, int i, int preNode, std::set<int> &path, std::vector< std::list<int> >&cycleLongest);
	int MergeCyclesToSuperNodes(vnl_matrix<unsigned char> &adj_matrix, std::vector< std::list<int> > &cycleList, std::vector< std::vector<std::list<int> > >&superNodeList, vnl_matrix<int> &superNodeConnection);
	void MergeCircles(std::list<int> &circle1, std::list<int> &circle2, vnl_vector<int> &commonNode, std::vector< std::list<int> > &leftCycleList);
	void BreakCircles(std::list<int> &circle, vnl_vector<int> &commonNode, std::vector< std::list<int> > &lines);
	template<class T> vnl_vector<T> VectorAnd(vnl_vector< T > const &v, vnl_vector< T > const &u);
	template<class T> vnl_vector<T> VectorOr(vnl_vector< T > const &v, vnl_vector< T > const &u);
	void GetConnectedLines(vnl_vector<int> &circle1, std::list<int> &circle2, vnl_vector<int> &commonNode, std::vector< std::list<int> > &lineList);
    int IsConnectedSuperNodeToNode(vnl_matrix<unsigned char> &adj_matrix, int node, std::list<int> &superNodePt);
    void CalculateCoordinatesForCircle(std::list<int> &circle, std::vector< std::list<int> > lineList, Point &center, double radius, std::vector<Point> &pointList);
	void GetLongestPath(std::vector< std::list<int>> &adjList, std::vector<int> &longestPath, std::vector<int> &tag, int node = 0);
	void depthSearch(std::vector< std::list<int>> &adjList, std::vector<int> &longestPath, std::vector<int> &tag, int node, std::vector<int> &bvisit);
	void depthAllSearch(std::vector< std::list<int>> &adjList, std::vector<int> &longestPath, std::vector<int> &tag, int node, std::vector<int> &bvisit);
	void GetAllNodeOnPath(std::vector< std::list<int>> &adjList, std::vector<int> &longestPath, std::vector<int> &tag, int node);
	void GetBackBones(std::vector< std::list<int>> &adjList, std::vector<int> backbones, std::vector<int> &tag, std::vector< std::pair<int, std::vector<int> > > &backboneMap);	
	void GetBackBonesMain(std::vector< std::list<int>> &adjList, std::vector<int> &backbones, std::vector< std::pair<int, std::vector<int> > > &backboneMap);

	protected slots:
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void HandleKeyPress(vtkObject* caller, long unsigned eventId, void* clientData, void* callData );
	void UpdateGraphView();
	
signals:
	void selection_Changed();

public:
	vtkSmartPointer<vtkGraphLayoutView> view;

private:
	vtkSmartPointer<vtkTable> dataTable;
	vnl_vector<double> colorVector;
	vnl_vector<double> featureColorVector;
	ObjectSelection *selection;
	ObjectSelection *selection2;
	   
	std::set<long int> colSelectIDs;
	
	QVTKWidget mainQTRenderWidget;
	//vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkTableToGraph> TTG;	
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkLookupTable> lookupTable;
	vtkSmartPointer<vtkLookupTable> edgeLookupTable;
	vtkSmartPointer<vtkCornerAnnotation> cornAnnotation;
	unsigned long observerTag;
	vnl_vector<double> edgeWeights;
	vnl_matrix<double> vertextList;

	std::map<int, int> indMapFromVertexToInd;
	std::vector<int> indMapFromIndToVertex;
	std::map<int, int> indMapFromVertexToClusInd;
	std::vector< std::vector<int> > indMapFromClusIndToVertex;
	std::vector< std::vector<int> > indMapFromClusIndToInd;
	std::map< std::pair< int, int>, int> edgeMapping;
	QString fileName;

	vnl_matrix<int> shortest_hop;
	std::vector< int> backbones;
	std::vector< std::pair<int, std::vector<int> > > chainList;

	// User Define progression
	int progressionStartID;
	int progressionEndID;
	bool bTrendStart;
	std::vector< int> progressionPath;
};

#endif
