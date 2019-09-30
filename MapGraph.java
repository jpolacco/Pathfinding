/**
 * @author UCSD MOOC development team and Joe Polacco
 *
 * Note to Airbnb team: take a look at my dijkstra and aStarSearch method implementations
 * starting at line 258.
 *
 * A class which represents a graph of geographic locations
 * Nodes in the graph are intersections between
 *
 * Class: MapGraph
 *
 * Modifications made to MapGraph (what and why):
 *
 *  I first wrote a confusing data structure called edgeLengths:
 *  private HashMap<GeographicPoint, HashMap<GeographicPoint, Double>>  edgeLengths;
 *
 *  I was finding this approach intractable, so I instead designed a MapEdge class and a MapNode class.
 *  This avoided having to use a HashMap of HashMap's.
 *
 *  My MapNode class wrapped both the GeographicPoint of the node and the edges from each GeographicPoint,
 *  a List of MapEdge's.
 *
 *  My MapEdge class encapsulated both the starting and ending GeographicPoint, as well
 *  as the road name and distance.
 *
 *  I created an instance variable called map that maps a GeographicPoint to a MapNode class.
 *  private Map<GeographicPoint, MapNode> map;
 *
 *  Thus, I no longer needed to map a GeographicPoint to a HashMap of GeometricPoint/Double key/value pairs, as the
 *  MapNode class kept track of the edges.  This more object oriented approach made the syntax a lot easier
 *  to read and made my code more readable.
 *
 *  Joe Polacco
 *  June, 2019
 *
 */
package roadgraph;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.function.Consumer;
import java.util.PriorityQueue;

import geography.GeographicPoint;
//import geography.RoadSegment;
import util.GraphLoader;



public class MapGraph {

	private Map<GeographicPoint, MapNode> map;
	private int numIntersections;
	private int numEdges;

	// avoided this ugly data structure by designing MapEdge and MapNode classes:
	// private HashMap<GeographicPoint, HashMap<GeographicPoint, Double>>  edgeLengths;

	/**
	 * Create a new empty MapGraph
	 */
	public MapGraph()
	{
		map = new HashMap<GeographicPoint, MapNode>();
		numIntersections = 0;
		numEdges = 0;
	}


	/**
	 * Get the number of vertices (road intersections) in the graph
	 * @return The number of vertices in the graph.
	 */
	public int getNumVertices()
	{
		return numIntersections;
	}

	/**
	 * Return the intersections, which are the vertices in this graph.
	 * @return The vertices in this graph as GeographicPoints
	 */
	public Set<GeographicPoint> getVertices()
	{
		return map.keySet();
	}

	/**
	 * Get the number of road segments in the graph
	 * @return The number of edges in the graph.
	 */
	public int getNumEdges()
	{
		return numEdges;
	}

	/** Add a node corresponding to an intersection at a Geographic Point
	 * If the location is already in the graph or null, this method does
	 * not change the graph.
	 * @param location  The location of the intersection
	 * @return true if a node was added, false if it was not (the node
	 * was already in the graph, or the parameter is null).
	 */
	public boolean addVertex(GeographicPoint location)
	{
		/*
		1.  if GeographicPoint already exists in the map or is null, then return false
		2.  else, we add the GeographicPoint to the map and then return true
		We have an instance variable called map.
		map maps GeographicPoints to a mapNode object.
		*/

		if (location == null || map.containsKey(location))
			return false;
		//create an intersection with an empyt list of edges
		MapNode node = new MapNode(location, new ArrayList<MapEdge>());
		map.put(location, node);
		numIntersections++;
		return true;
	}

	/**
	 * Adds a directed edge to the graph from pt1 to pt2.
	 * Precondition: Both GeographicPoints have already been added to the graph
	 * @param from The starting point of the edge
	 * @param to The ending point of the edge
	 * @param roadName The name of the road
	 * @param roadType The type of the road
	 * @param length The length of the road, in km
	 * @throws IllegalArgumentException If the points have not already been
	 *   added as nodes to the graph, if any of the arguments is null,
	 *   or if the length is less than 0.
	 */
	public void addEdge(GeographicPoint from, GeographicPoint to, String roadName,
						String roadType, double length) throws IllegalArgumentException {
		//TODO: Implement this method
		if (  !map.containsKey(from) || !map.containsKey(to) || roadName == null
				|| roadType == null || length < 0.0)
			throw new IllegalArgumentException();
		MapEdge edge = new MapEdge(from, to, roadName, length);
		MapNode node = map.get(from);
		node.addEdge(edge);
		numEdges++;

	}

	/** Find the path from start to goal using breadth first search
	 *
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		Consumer<GeographicPoint> temp = (x) -> {};
		return bfs(start, goal, temp);
	}

	/** Find the path from start to goal using breadth first search
	 *
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start,
									 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		MapNode startNode = map.get(start);
		MapNode endNode = map.get(goal);
		if (startNode == null || endNode == null) {
			System.out.println("Start or goal node is null!  No path exists.");
			return new LinkedList<GeographicPoint>();
		}

		HashSet<GeographicPoint> visited = new HashSet<GeographicPoint>();
		Queue<MapNode> toExplore = new LinkedList<MapNode>();
		HashMap<GeographicPoint, GeographicPoint> parentMap = new HashMap<GeographicPoint, GeographicPoint>();

		toExplore.add(startNode);
		boolean found = false;
		while (!toExplore.isEmpty()){
			MapNode curr = toExplore.remove();

			//used for visualization
			nodeSearched.accept(curr.getGeographicPoint());

			if (curr == endNode){ // maybe need .equals here?
				found = true;
				break;
			}
			List<GeographicPoint> neighbors = curr.getNeighbors();
			for (GeographicPoint neighbor : neighbors){
				if (!visited.contains(neighbor)){
					visited.add(neighbor);
					parentMap.put(neighbor, curr.getGeographicPoint());
					toExplore.add(map.get(neighbor));

				}
			}
		}

		if (!found){
			System.out.println("No path exists");
			//return new ArrayList<GeographicPoint>();
			return null;
		}

		//reconstruct the path
		LinkedList<GeographicPoint> path = new LinkedList<GeographicPoint>();
		GeographicPoint curr = goal;
		while (!curr.toString().equals(start.toString()) ){
			//for (int i = 0; i < 4; i++){
			path.addFirst(curr);
			curr = parentMap.get(curr);
		}
		path.addFirst(start);
		return path;

		// Hook for visualization.
		//nodeSearched.accept(next.getLocation());

		//nodeSearched.accept(next);
		//return null;
	}


	/** Finds the path from start to goal using Dijkstra's algorithm
	 *
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		Consumer<GeographicPoint> temp = (x) -> {};
		return dijkstra(start, goal, temp);
	}

	/** Finds the path from start to goal using Dijkstra's algorithm
	 *
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.
	 * @return The list of intersections that form the shortest path from
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start,
										  GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		//nodeSearched.accept(next.getLocation());
		MapNode startNode = map.get(start);
		MapNode endNode = map.get(goal);
		if (startNode == null || endNode == null) {
			System.out.println("Start or goal node is null!  No path exists.");
			return new LinkedList<GeographicPoint>();
		}

		HashSet<GeographicPoint> visited = new HashSet<GeographicPoint>();
		PriorityQueue<NodeDistanceFromStart> toExplore = new PriorityQueue<NodeDistanceFromStart>();
		HashMap<GeographicPoint, GeographicPoint> parentMap = new HashMap<GeographicPoint, GeographicPoint>();
		HashMap<GeographicPoint, Double> distanceFromStart = new HashMap<GeographicPoint, Double>();
		//build up Priority Queue
		NodeDistanceFromStart pqNode = new NodeDistanceFromStart(startNode, 0.0);
		toExplore.add(pqNode);
		//add to distanceFromStart map
		distanceFromStart.put(start, 0.0);

		//add all other nodes in the map as distances of infinity from start
		for (Map.Entry<GeographicPoint, MapNode> entry : map.entrySet()) {
			pqNode = new NodeDistanceFromStart(map.get(entry.getKey()), Double.POSITIVE_INFINITY);
			toExplore.add(pqNode);
		}

		boolean found = false;

		int nodesExplored = 0;  // used for debugging

		while (!toExplore.isEmpty()){

			NodeDistanceFromStart currNodeDistanceFromStart = toExplore.poll();
			//check to see if polled element in queue has already been visited
			if (visited.contains(currNodeDistanceFromStart.node.getGeographicPoint()))
				continue;
			nodesExplored++;

			MapNode curr = currNodeDistanceFromStart.node;

			//used for visualization
			nodeSearched.accept(curr.getGeographicPoint());

			if (curr == endNode){ // maybe need .equals here?
				found = true;
				break;
			}
			List<GeographicPoint> neighbors = curr.getNeighbors();
			for (GeographicPoint neighbor : neighbors){
				if (!visited.contains(neighbor)){
					//update new distance from start for each neighbor not in visited (only if needed)
					double distance = currNodeDistanceFromStart.distanceFromStart + neighbor.distance(curr.getGeographicPoint());

					// 1. Add neighbor to priority queue
					// okay if adding duplicating nodes -- priority queue will poll the lowest distance from start
					toExplore.add(new NodeDistanceFromStart(map.get(neighbor),distance));


					Double currentDistance = distanceFromStart.get(neighbor);
					if (currentDistance == null || distance < currentDistance ){

						//2. Update hash map with shortest distances from start:
						//   add new entry into map or override old entry
						distanceFromStart.put(neighbor,distance);
						// 3. update parent map
						parentMap.put(neighbor,curr.getGeographicPoint());

					}
					// else do nothing -- old distance is still shorter or equal

				} // end if not visited
			}
			//4. Add current Geographic Point to visited
			visited.add(curr.getGeographicPoint());
		}

		if (!found){
			System.out.println("No path exists");
			//return new ArrayList<GeographicPoint>();
			return null;
		}

		//reconstruct the path
		LinkedList<GeographicPoint> path = new LinkedList<GeographicPoint>();
		GeographicPoint curr = goal;

		System.out.println("Dijkstra Search explored " + nodesExplored + " nodes");

		while (!curr.toString().equals(start.toString()) ){
			//for (int i = 0; i < 4; i++){
			path.addFirst(curr);
			curr = parentMap.get(curr);
		}
		path.addFirst(start);
		return path;

		//return null;
	}

	/** Find the path from start to goal using A-Star search
	 *
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		Consumer<GeographicPoint> temp = (x) -> {};
		return aStarSearch(start, goal, temp);
	}

	/** Find the path from start to goal using A-Star search
	 *
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.
	 * @return The list of intersections that form the shortest path from
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start,
											 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		// Hook for visualization.
		//nodeSearched.accept(next.getLocation());
		MapNode startNode = map.get(start);
		MapNode endNode = map.get(goal);
		if (startNode == null || endNode == null) {
			System.out.println("Start or goal node is null!  No path exists.");
			return new LinkedList<GeographicPoint>();
		}

		HashSet<GeographicPoint> visited = new HashSet<GeographicPoint>();
		PriorityQueue<NodeDistanceFromStart> toExplore = new PriorityQueue<NodeDistanceFromStart>();
		HashMap<GeographicPoint, GeographicPoint> parentMap = new HashMap<GeographicPoint, GeographicPoint>();
		HashMap<GeographicPoint, Double> distanceFromStart = new HashMap<GeographicPoint, Double>();
		//build up Priority Queue
		NodeDistanceFromStart pqNode = new NodeDistanceFromStart(startNode, 0.0);
		toExplore.add(pqNode);
		//add to distanceFromStart map
		distanceFromStart.put(start, 0.0);
		//add all other nodes in the map as distances of infinity from start
		for (Map.Entry<GeographicPoint, MapNode> entry : map.entrySet()) {
			pqNode = new NodeDistanceFromStart(map.get(entry.getKey()), Double.POSITIVE_INFINITY);
			toExplore.add(pqNode);
		}

		boolean found = false;

		int nodesExplored = 0;  // used for debugging
		while (!toExplore.isEmpty()){

			NodeDistanceFromStart currNodeDistanceFromStart = toExplore.poll();

			//check to see if polled element in queue has already been visited
			if (visited.contains(currNodeDistanceFromStart.node.getGeographicPoint()))
				continue;
			nodesExplored++;
			MapNode curr = currNodeDistanceFromStart.node;

			//used for visualization
			nodeSearched.accept(curr.getGeographicPoint());

			if (curr == endNode){ // maybe need .equals here?
				found = true;
				break;
			}
			List<GeographicPoint> neighbors = curr.getNeighbors();
			for (GeographicPoint neighbor : neighbors){
				if (!visited.contains(neighbor)){
					//update new distance from start for each neighbor not in visited (only if needed)
					double distance = currNodeDistanceFromStart.distanceFromStart + neighbor.distance(curr.getGeographicPoint());

					// 1. Add neighbor to priority queue
					// okay if adding duplicating nodes -- priority queue will poll the lowest distance from start
					//A* heuristic: calculate Euclidean distance from current node to goal node
					double heuristicDistance = neighbor.distance(goal);

					toExplore.add(new NodeDistanceFromStart(map.get(neighbor),distance + heuristicDistance));


					Double currentDistance = distanceFromStart.get(neighbor);
					if (currentDistance == null || distance < currentDistance ){

						//2. Update hash map with shortest distances from start:
						//   add new entry into map or override old entry
						distanceFromStart.put(neighbor,distance);
						// 3. update parent map
						parentMap.put(neighbor,curr.getGeographicPoint());

					}
					// else do nothing -- old distance is still shorter or equal

				} // end if not visited
			}
			//4. Add current Geographic Point to visited
			visited.add(curr.getGeographicPoint());
		}

		if (!found){
			System.out.println("No path exists");
			//return new ArrayList<GeographicPoint>();
			return null;
		}

		System.out.println("A* Search explored " + nodesExplored + " nodes");
		//reconstruct the path
		LinkedList<GeographicPoint> path = new LinkedList<GeographicPoint>();
		GeographicPoint curr = goal;
		while (!curr.toString().equals(start.toString()) ){
			//for (int i = 0; i < 4; i++){
			path.addFirst(curr);
			curr = parentMap.get(curr);
		}
		path.addFirst(start);
		return path;
		//return null;
	}

	/** Solves the Traveling Salesman Problem using a Greedy Algorithm
	 *
	 * @param start The starting AND ending location
	 * @param destinations All the intersections we need to visit exactly once
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the ordered path that visits all intersections from start and ending back at start
	 */
	public List<GeographicPoint> travelingSalesman(GeographicPoint start, List<GeographicPoint> destinations, Consumer<GeographicPoint> nodeSearched) {


		LinkedList<GeographicPoint> shortestPathWithIntermediateIntersections = new LinkedList<>();
		double shortestDistance = Double.MAX_VALUE;
		int sizeOfDestinations = destinations.size();

		// Run traveling salesman n times, where n is the number of destinations
		// On the ith run of TSP, choose a random edge for the ith node instead of the  shortest distance
		// and keep track of the shortest path to see if this heuristic improves upon greedy
		for (int i = 0; i < destinations.size(); i++) {
			//clone destinations to reuse it for the ith pass
			List<GeographicPoint> destinationsClone = new ArrayList<>();
			for (GeographicPoint gp : destinations )
				destinationsClone.add(gp);

			HashSet<GeographicPoint> visited = new HashSet<>();
			//List<GeographicPoint> destination = destinations.clone();
			LinkedList<GeographicPoint> path = new LinkedList<>();
			LinkedList<GeographicPoint> pathWithIntermediateIntersections = new LinkedList<>();

			path.add(start);  // add the first intersection.  It is our starting point
			pathWithIntermediateIntersections.add(start);
			GeographicPoint current = start;

			double totalPathDistance = 0.0;
			while (destinations.size() > 0) {

				//find closest intersection in route that has not yet been visited(greedy)
				List<GeographicPoint> closestPathOption = new LinkedList<>(); // used to keep track of intermediate intersections
				double closestIntersectionDistance = Double.POSITIVE_INFINITY;
				GeographicPoint closestIntersection = destinations.get(0);
				if ((int) (Math.random()*sizeOfDestinations) == 1) {
					//System.out.println("hi");
					List<GeographicPoint> pathOption = new LinkedList<GeographicPoint>();
					pathOption = aStarSearch(current, destinations.get(0));
					double distance = distanceAlongPath(pathOption);
					closestIntersectionDistance = distance;
					closestIntersection = destinations.get(0);
					closestPathOption = pathOption;
				}
				else {

					/*
					if (distance < closestIntersectionDistance) {
						closestIntersectionDistance = distance;
						closestIntersection = destinations.get(j);
						closestPathOption = pathOption;
					}
					*/

					for (int j = 0; j < destinations.size(); j++) {


						List<GeographicPoint> pathOption = new LinkedList<GeographicPoint>();
						pathOption = aStarSearch(current, destinations.get(j));
						double distance = distanceAlongPath(pathOption);
						if (distance < closestIntersectionDistance) {
							closestIntersectionDistance = distance;
							closestIntersection = destinations.get(j);
							closestPathOption = pathOption;
						}
					}
				}
				/*
				for (GeographicPoint gp : destinations) {

					List<GeographicPoint> pathOption = new LinkedList<GeographicPoint>();
					pathOption = aStarSearch(current, gp);
					double distance = distanceAlongPath(pathOption);
					if (distance < closestIntersectionDistance) {
						closestIntersectionDistance = distance;
						closestIntersection = gp;
						closestPathOption = pathOption;
					}
				}
				*/

				totalPathDistance += closestIntersectionDistance;
				path.add(closestIntersection);
				closestPathOption.remove(0);
				pathWithIntermediateIntersections.addAll(closestPathOption);
				destinations.remove(closestIntersection);
				current = closestIntersection;

			}
			path.add(start);  // complete the loop

			System.out.println(pathWithIntermediateIntersections);
			System.out.println("Total distance of this path: " + totalPathDistance);

			if (totalPathDistance < shortestDistance) {
				shortestDistance = totalPathDistance;
				shortestPathWithIntermediateIntersections = pathWithIntermediateIntersections;
			}

			destinations = destinationsClone;
			//return path;
			//return pathWithIntermediateIntersections;
		} //end outer for
		return shortestPathWithIntermediateIntersections;
	}

	/**
	 *
	 * @param path The list of intersections that make up a path
	 * @return  Return the total distance travelled to visit all intersections in the path
	 */
	private double distanceAlongPath(List<GeographicPoint> path){
		double distance = 0.0;
		GeographicPoint current = path.get(0);
		for (int i = 1; i < path.size(); i++){
			distance += current.distance(path.get(i));
			current = path.get(i);
		}
		return distance;
	}

	public static void main(String[] args)
	{
		System.out.print("Making a new map...");
		MapGraph firstMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", firstMap);
		System.out.println("DONE.");

		/*

		System.out.println("Testing Traveling Salesman Problem");
		List<GeographicPoint> nodes = new ArrayList<>();
		nodes.add(new GeographicPoint(4,1));
		nodes.add(new GeographicPoint(4,2));
		nodes.add(new GeographicPoint(5,1));
		nodes.add(new GeographicPoint(7,3));
		nodes.add(new GeographicPoint(8,-1));
		nodes.add(new GeographicPoint(7,3));
		nodes.add(new GeographicPoint(4,0));
		System.out.println(firstMap.travelingSalesman(new GeographicPoint(1,1), nodes,null));
		*/


		System.out.println("BFS testing");
		//GeographicPoint testStart = new GeographicPoint(1.0, 1.0);
		//GeographicPoint testEnd = new GeographicPoint(8.0, -1.0);
		GeographicPoint testStart = new GeographicPoint(4.0, -1.0);
		GeographicPoint testEnd = new GeographicPoint(4.0, -1.0);
		System.out.print(firstMap.bfs(testStart, testEnd));
		MapGraph testMap = new MapGraph();
		GraphLoader.loadRoadMap("data/maps/utc.map", testMap);

		// A very simple test using real data
		testStart = new GeographicPoint(32.869423, -117.220917);
		testEnd = new GeographicPoint(32.869255, -117.216927);
		System.out.println("Test 2 using bfs: ");
		System.out.println(testMap.bfs(testStart,testEnd));



		MapGraph simpleTestMap = new MapGraph();
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", simpleTestMap);

		testStart = new GeographicPoint(1.0, 1.0);
		testEnd = new GeographicPoint(8.0, -1.0);

		System.out.println("Test 1 using simpletest: Dijkstra should be 9 and AStar should be 5");
		List<GeographicPoint> testroute = simpleTestMap.dijkstra(testStart,testEnd);
		List<GeographicPoint> testroute2 = simpleTestMap.aStarSearch(testStart,testEnd);


		testMap = new MapGraph();
		GraphLoader.loadRoadMap("data/maps/utc.map", testMap);

		// A very simple test using real data
		testStart = new GeographicPoint(32.869423, -117.220917);
		testEnd = new GeographicPoint(32.869255, -117.216927);
		System.out.println("Test 2 using utc: Dijkstra should be 13 and AStar should be 5");
		testroute = testMap.dijkstra(testStart,testEnd);
		testroute2 = testMap.aStarSearch(testStart,testEnd);


		// A slightly more complex test using real data
		testStart = new GeographicPoint(32.8674388, -117.2190213);
		testEnd = new GeographicPoint(32.8697828, -117.2244506);
		System.out.println("Test 3 using utc: Dijkstra should be 37 and AStar should be 10");
		testroute = testMap.dijkstra(testStart,testEnd);
		testroute2 = testMap.aStarSearch(testStart,testEnd);
		System.out.println(testroute2);
		System.out.println(testMap.distanceAlongPath(testroute2));



		MapGraph theMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/maps/utc.map", theMap);
		System.out.println("DONE.");

		GeographicPoint start = new GeographicPoint(32.8648772, -117.2254046);
		GeographicPoint end = new GeographicPoint(32.8660691, -117.217393);


		List<GeographicPoint> route = theMap.dijkstra(start,end);
		List<GeographicPoint> route2 = theMap.aStarSearch(start,end);



	}

}

// class used for Priority Queue in Dijkstra's Algorithm
class NodeDistanceFromStart implements Comparable<NodeDistanceFromStart>{
	//public String word;
	//
	// code from Dijkstra:
	//PriorityQueue<NodeDistanceFromStart> toExplore = new PriorityQueue<NodeDistanceFromStart>();
	public GeographicPoint gp;
	public MapNode node;
	public double distanceFromStart;

	public NodeDistanceFromStart(MapNode m, double n){
		node = m;
		distanceFromStart = n;
	}

	public int compareTo(NodeDistanceFromStart that){
		Double distance = this.distanceFromStart - that.distanceFromStart;
		if (distance > 0)
			return 1;
		else if (distance < 0)
			return -1;
		else
			return 0;
	}

	public String toString(){
		return node + " " + distanceFromStart;
	}
}

