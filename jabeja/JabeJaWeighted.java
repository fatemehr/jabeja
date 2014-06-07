package edu.cmu.graphchi.apps.jabeja;

import edu.cmu.graphchi.*;
import edu.cmu.graphchi.datablocks.IntArrayConverter;
import edu.cmu.graphchi.engine.GraphChiEngine;
import edu.cmu.graphchi.engine.VertexInterval;
import edu.cmu.graphchi.preprocessing.EdgeProcessor;
import edu.cmu.graphchi.preprocessing.FastSharder;
import edu.cmu.graphchi.preprocessing.VertexProcessor;


import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

/**
 * implements JabeJa for weighted graphs
 * Iteratively computes a balanced partitioning of a given graph
 * by performing local operations at each vertex.
 * @author Fatemeh Rahimian
 */
public class JabeJaWeighted implements GraphChiProgram<Integer[], Integer[]> {
    private static Random rnd = new Random(456);
    private static Logger logger = ChiLogger.getLogger("jabeja");
    private static int numPartitions;
    private static float TEMPERATURE = (float) 1.5;			// if set to 1, it means no simulated annealing is applied
    private static float TEMPERATUREDelta = (float) 0.003;
	
    private static enum MType {
		INFO, ACK, NACK, SWAP
	}
    private static final int ISWAITING = 0;
    private static final int SENDER = 0;
    private static final int COLOR = 1;
    private static final int TYPE = 2;
    private static int WEIGHT; 

    private static Integer[] colorInventory;
    
    private static AtomicInteger cut;
	private static MessageRelay mail = new MessageRelay();

    // -----------------------------------------------------------------
    public void update(ChiVertex<Integer[], Integer[]> vertex, GraphChiContext context)  {
        if (context.getIteration() == 0) {		/* Initialize on first iteration */
        	int vertexColor;
        	int isVertexWaiting;
			vertexColor = rnd.nextInt(numPartitions);
			while (colorInventory[vertexColor] > context.getNumVertices()/numPartitions)
				vertexColor = (vertexColor + 1) % numPartitions;
			colorInventory[vertexColor]++;
			
			isVertexWaiting = 0;
			Integer[] nodeValue = { isVertexWaiting, vertexColor };
            vertex.setValue(nodeValue);
            
			// format [sender, sender color, MType, color distribution]
            for(int i=0; i<vertex.numOutEdges(); i++) {
            	Integer[] info = vertex.getOutEdgeValue(i);
				info[SENDER] = vertex.getId();
				info[TYPE] = MType.INFO.ordinal();
				info[COLOR] = vertexColor;
				
                vertex.outEdge(i).setValue(info);
//                logger.info("iteration: " + context.getIteration() + " vertex:" + vertex.getId() + " color: " + vertex.getOutEdgeValue(i)[COLOR]+ " edge to neighbor i=" + i + " "+ vertex.getOutEdgeValue(i)[WEIGHT]);
            }
            
            
            
        } else {	// for all other iterations
        	int inDegree = vertex.numInEdges();
        	int outDegree = vertex.numOutEdges();
        	int weightedInDegree = 0;
        	
        	// node to inEdge/outEdge mapping <VertexID , edgeIndex>
        	HashMap<Integer, Integer> inNeighborsMap = new HashMap<Integer, Integer>();
        	HashMap<Integer, Integer> outNeighborsMap = new HashMap<Integer, Integer>();

        	HashMap<Integer, Integer[]> inEdges = new HashMap<Integer, Integer[]>();
        	HashMap<Integer, Integer[]> outEdges = new HashMap<Integer, Integer[]>();
			
        	// STEP 1: read all the in/out values
        	Integer[] nodeValue = vertex.getValue();
        	
        	int[] population = new int[numPartitions];
        	Arrays.fill(population, 0);
        	
        	for (int i = 0; i < inDegree; i++) {
        		inNeighborsMap.put(vertex.inEdge(i).getVertexId(), i);
				inEdges.put(i, vertex.inEdge(i).getValue());
        		population[vertex.inEdge(i).getValue()[COLOR]] += vertex.inEdge(i).getValue()[WEIGHT];
        		weightedInDegree += vertex.inEdge(i).getValue()[WEIGHT];
			}
        	for (int i = 0; i < outDegree; i++) {
        		outNeighborsMap.put(vertex.getOutEdgeId(i), i);
				outEdges.put(i, vertex.getOutEdgeValue(i));
			}
        	
        	
        	// read from mail (this does not update the population)
        	ArrayList<Integer[]> messages = new ArrayList<Integer[]>(mail.get(vertex.getId()));
        	
        	
        	// STEP 2: process / update state
        	ArrayList<Integer[]> swapRequests = new ArrayList<Integer[]>();
        	Integer[] inValue;
        	Integer[] outValue;
        	int outIndex;
        	
        	for (Integer[] mailreq : messages) {
//        		logger.info("received a mail from node:" + mailreq[SENDER] + " with type:" + mailreq[TYPE]);
        		MType mType = MType.values()[mailreq[TYPE]];
        		if (mType == MType.ACK) {
        			nodeValue[COLOR] = mailreq[COLOR];
        			nodeValue[ISWAITING] = 0;			// value 0 indicates that isVertexWaiting = false;
        			vertex.setValue(nodeValue);
        		} else if (mType == MType.NACK) {
        			nodeValue[ISWAITING] = 0;			// value 0 indicates that isVertexWaiting = false;
        			vertex.setValue(nodeValue);
        		}
        		else if (mType == MType.SWAP) {
        			swapRequests.add(mailreq);
        		}
        	}
        	
        	for (int i = 0; i < inDegree; i++) {
        		inValue = inEdges.get(i);
        		MType inType = MType.values()[inValue[TYPE]];

        		outIndex = outNeighborsMap.get(inValue[SENDER]);

        		if (inType == MType.INFO) {
        			// no state update, just reply with INFO
        			outValue = outEdges.get(outIndex);
        			outValue[TYPE] = MType.INFO.ordinal();
        			outEdges.put(outIndex, outValue);
        		}
        		else if (inType == MType.ACK) {
        			// change your color, but before that make corrections to population array
        			population[nodeValue[COLOR]] += inValue[WEIGHT];
        			population[inValue[COLOR]] -= inValue[WEIGHT];
        			nodeValue[COLOR] = inValue[COLOR];
        			nodeValue[ISWAITING] = 0;	// value 0 indicates that isVertexWaiting = false;
                    vertex.setValue(nodeValue);
                    
        			outValue = outEdges.get(outIndex);
        			outValue[TYPE] = MType.INFO.ordinal();
        			outEdges.put(outIndex, outValue);
        		}
        		else if (inType == MType.NACK) {
        			nodeValue[ISWAITING] = 0;
                    vertex.setValue(nodeValue);
                    
        			outValue = outEdges.get(outIndex);
        			outValue[TYPE] = MType.INFO.ordinal();
        			outEdges.put(outIndex, outValue);
        		}
        		else if (inType == MType.SWAP) {
        			if (nodeValue[ISWAITING] == 0) {
        				swapRequests.add(inValue);
        			}
			        else {
	        			outValue = outEdges.get(outIndex);
	        			outValue[TYPE] = MType.NACK.ordinal();
	        			outEdges.put(outIndex, outValue);
			        }
        		}
        		
        	}
        	
        	boolean colorChanged = false;
        	// process the swap requests, if any
    		if (nodeValue[ISWAITING] == 0 && swapRequests.size() > 0) {
//    			logger.info(vertex.getId() + " is proccessing " + swapRequests.size() + " swap request(s).");
    			Collections.shuffle(swapRequests, rnd);			// shuffle all the incoming requests
				for (Integer[] req: swapRequests) {
					if (inNeighborsMap.containsKey(req[SENDER])) { // the request is sent by a neighbor, so update the edge value
						if (colorChanged == false && swap(vertex.getValue()[COLOR], population, req, true) > 0) {
							outIndex = outNeighborsMap.get(req[SENDER]);
		                    Integer[] edgeValue = vertex.outEdge(outIndex).getValue();
		        			edgeValue[TYPE] = MType.ACK.ordinal();
		        			edgeValue[COLOR] = nodeValue[COLOR];
		        			outEdges.put(outIndex, edgeValue);
		        			
		        			nodeValue[COLOR] = req[COLOR];
		        			nodeValue[ISWAITING] = 0;
		        			vertex.setValue(nodeValue);
		        			
//		        			logger.info("****** Swap request is accepted");
		        			
		        			colorChanged = true;
						} 
						else {
							outIndex = outNeighborsMap.get(req[SENDER]);
	        				Integer[] edgeValue = vertex.outEdge(outIndex).getValue();
		        			edgeValue[TYPE] = MType.NACK.ordinal();
		        			outEdges.put(outIndex, edgeValue);

						}
					}
					else {	// the request is sent by a random node
						if (colorChanged == false && swap(vertex.getValue()[COLOR], population, req, false) > 0) {
	        				Integer[] response = new Integer[4 + numPartitions];
	        				Arrays.fill(response, -1);
	        				response[SENDER] = vertex.getId();
	        				response[TYPE] = MType.ACK.ordinal();	 
	        				response[COLOR] = nodeValue[COLOR];
	        				mail.send(req[SENDER], response);
//	        				logger.info("sending an ACK response to node " + req[SENDER]);
	        				
	        				nodeValue[COLOR] = req[COLOR];
		        			nodeValue[ISWAITING] = 0;
		        			vertex.setValue(nodeValue);
		        			
		        			colorChanged = true;
						}
						else {
	        				Integer[] response = new Integer[4 + numPartitions];
	        				Arrays.fill(response, -1);
	        				response[SENDER] = vertex.getId();
	        				response[TYPE] = MType.NACK.ordinal();
	        				mail.send(req[SENDER], response);
//	        				logger.info("sending a NACK response to node " + req[SENDER]);

						}
					}
						
				}
    		}
    		else if (swapRequests.size() > 0) {
//    			logger.info(vertex.getId() + " has to reject " + swapRequests.size() + " swap request(s).");
				for (Integer[] req: swapRequests) {
					if (inNeighborsMap.containsKey(req[SENDER])) { // the request is sent by a neighbor, so update the edge value
						outIndex = outNeighborsMap.get(req[SENDER]);
	    				Integer[] edgeValue = vertex.outEdge(outIndex).getValue();
	        			edgeValue[TYPE] = MType.NACK.ordinal();
	        			outEdges.put(outIndex, edgeValue);
					}
					else {
        				Integer[] response = new Integer[4 + numPartitions];
        				Arrays.fill(response, -1);
        				response[SENDER] = vertex.getId();
        				response[TYPE] = MType.NACK.ordinal();
        				mail.send(req[SENDER], response);
//        				logger.info("sending a NACK response to node " + req[SENDER]);
					}
				}
    		}
        	
        	// Send a swap request, either to a neighbor or to a random node in the graph
    		// first, try to find the best neighbor to swap with, if no good local swap is possible select a random node
        	if (nodeValue[ISWAITING] == 0) {
        		int bestNeighbor;
        		if (outDegree > 0 && (bestNeighbor = selectBestNeighbor(vertex.getValue()[COLOR], population, outEdges)) != -1 ) {
	        		Integer[] edgeValue = vertex.getOutEdgeValue(bestNeighbor);
	        		edgeValue[SENDER] = vertex.getId();
	        		edgeValue[COLOR] = nodeValue[COLOR];
	        		edgeValue[TYPE] = MType.SWAP.ordinal();
	        		for (int i = 0; i < numPartitions; i++)
	        			edgeValue[3 + i] = population[i];
	        		
//	        		logger.info("sending swap request... to neighbor:" + bestNeighbor);
	        		outEdges.put(bestNeighbor, edgeValue);
	        		nodeValue[ISWAITING] = 1;

        		}
        		else {
        			Integer[] req = new Integer[4 + numPartitions];
	        		req[SENDER] = vertex.getId();
	        		req[COLOR] = nodeValue[COLOR];
	        		req[TYPE] = MType.SWAP.ordinal();
	        		for (int i = 0; i < numPartitions; i++)
	        			req[3 + i] = population[i];
	        		
	        		Long N = context.getNumVertices();
	        		Integer randomNodeId = rnd.nextInt(N.intValue());
	        		if (vertex.getId() != randomNodeId && !inNeighborsMap.containsKey(randomNodeId)) {
	        			mail.send(randomNodeId, req);
		        		nodeValue[ISWAITING] = 1;
	        		}
//	        		logger.info("sending swap request... to random node id:" + randomNodeId);
        		}
        		vertex.setValue(nodeValue);
        	}
        	
        	
        	// write the corresponding out-values
        	for (int i = 0; i < outDegree; i++) {
        		Integer[] edgeValue = outEdges.get(i);
        		
        		for (int j= 0; j < numPartitions; j++)
        			edgeValue[3 + j] = population[j];
        		if (edgeValue[TYPE] != MType.ACK.ordinal()) // if you are not sending an ACK send your current color to the neighbor, in case of ACK let it be the old color
        			edgeValue[COLOR] = nodeValue[COLOR];
        		
       			vertex.outEdge(i).setValue(edgeValue);
//              logger.info("iteration: " + context.getIteration() + " vertex:" + vertex.getId() + " color: " + vertex.getOutEdgeValue(i)[COLOR]+ " edge to neighbor i=" + i + " "+ vertex.getOutEdgeValue(i)[WEIGHT]);
//        		logger.info("outedgeValue for neighbor:" + i + " is [sender:" + edgeValue[SENDER] + " ,type:" + edgeValue[TYPE] + "]");
        	}
        	
        	
//        	logger.info("node id: " + vertex.getId() + " round: " + context.getIteration() + " color: " + vertex.getValue()[COLOR]);
        	cut.addAndGet(weightedInDegree - population[vertex.getValue()[COLOR]]);
        }
    }

	// -----------------------------------------------------------------
    private int selectBestNeighbor(int selfColor, int[] localPopulation, HashMap<Integer, Integer[]> outEdges) {
    	int bestId = -1;
    	double best = 0;
    	double benefit;
    	
		for (int id : outEdges.keySet()) {
			benefit = swap(selfColor, localPopulation, outEdges.get(id), true);
			if (benefit > best) {
				best = benefit;
				bestId = id;
			}
		}

		return bestId;
	}

	// -----------------------------------------------------------------
	private double swap(int selfColor, int[] localPopulation, Integer[] inValue, boolean isNeighbor) {
		float coefficient;
		
		int tempNodeColor = inValue[COLOR];
		
		if (tempNodeColor == selfColor)
			return 0;
		
		int c1SelfNodeBenefit = localPopulation[selfColor];
		int c1TempNodeBenefit = inValue[3+tempNodeColor];	// the population of current color
	
		int c2SelfNodeBenefit = localPopulation[tempNodeColor];
		int c2TempNodeBenefit = inValue[3+selfColor];		// the population of (possibly) future color
		
		if (isNeighbor == true) {	// if the swap of two neighbors are considered, the expected populations are overestimated by one for each node
			c2SelfNodeBenefit--;
			c2TempNodeBenefit--;
			coefficient = 1;
		}
		else
			coefficient = TEMPERATURE;
		
		double oldBenefit = Math.pow(c1SelfNodeBenefit, 2) + Math.pow(c1TempNodeBenefit, 2);
		double newBenefit = Math.pow(c2SelfNodeBenefit, 2) + Math.pow(c2TempNodeBenefit, 2);
		
		
		return newBenefit * coefficient - oldBenefit;
	}
	
    // -----------------------------------------------------------------
    public void beginIteration(GraphChiContext ctx) {
    	TEMPERATURE -= TEMPERATUREDelta;
    	if (TEMPERATURE < 1)
    		TEMPERATURE = 1;
    	
    	cut = new AtomicInteger(0);
    }

    public void endIteration(GraphChiContext ctx) {
        logger.setLevel(Level.INFO);
    	logger.info("iteration: " + ctx.getIteration() + " edge cut = " + cut.intValue() / 2);
//        logger.setLevel(Level.OFF);

    }

    public void beginInterval(GraphChiContext ctx, VertexInterval interval) {}

    public void endInterval(GraphChiContext ctx, VertexInterval interval) {}

    public void beginSubInterval(GraphChiContext ctx, VertexInterval interval) {}

    public void endSubInterval(GraphChiContext ctx, VertexInterval interval) {}

    /**
     * Initialize the sharder-program.
     * @param graphName
     * @param numShards
     * @return
     * @throws java.io.IOException
     */
    protected static FastSharder<Integer[], Integer[]> createSharder(String graphName, int numShards) throws IOException {
        return new FastSharder<Integer[], Integer[]>(graphName, numShards, new VertexProcessor<Integer[]>() {
            public Integer[] receiveVertexValue(int vertexId, String token) {
                return new Integer[2];
            }
        }, new EdgeProcessor<Integer[]>() {
            public Integer[] receiveEdge(int from, int to, String token) {
            	Integer[] val = new Integer[4 + numPartitions];
            	for (int i = 0; i < 3 + numPartitions; i++)
            		val[i] = 0;
            	val[3 + numPartitions] = Integer.parseInt(token);
                return val;
            }
        }, new IntArrayConverter(2), new IntArrayConverter(numPartitions + 4));
    }


    /**
     * Usage: java edu.cmu.graphchi.demo.ConnectedComponents graph-name num-shards filetype(edgelist|adjlist)
     * For specifying the number of shards, 20-50 million edges/shard is often a good configuration.
     */
    public static void main(String[] args) throws  Exception {
    	String fileType = "edgelist";
    	String baseFilename = "";
    	int nShards = 1;
    	int numberOfRounds = 10;
    	try {
	        baseFilename = args[0];
	        nShards = Integer.parseInt(args[1]);
	        numberOfRounds = Integer.parseInt(args[2]);
	        numPartitions = Integer.parseInt(args[3]);
	        TEMPERATURE = (args.length > 4 ? Float.parseFloat(args[4]) : (float) 1.5);
	        
	        WEIGHT= 3 + numPartitions;
	        colorInventory = new Integer[numPartitions];
    	}
    	catch (Exception e) {
    		logger.info(e.toString());
    		logger.warning("The input arguments are not properly defined.");
    		logger.info("... edu.cmu.graphchi.apps.JabeJaWeighted <input-file> <#shards> <#rounds> <#partitions> <[optional: initial temperature]>");
    		System.exit(0);
    	}
    	
        logger.addHandler(new FileHandler(baseFilename + "_log.txt"));
        
        /* Create shards */
        FastSharder sharder = createSharder(baseFilename, nShards);
        if (baseFilename.equals("pipein")) {     // Allow piping graph in
            sharder.shard(System.in, fileType);
        } else {
            if (!new File(ChiFilenames.getFilenameIntervals(baseFilename, nShards)).exists()) {
                sharder.shard(new FileInputStream(new File(baseFilename)), fileType);
            } else {
                logger.info("Found shards -- no need to preprocess");
            }
        }

        Arrays.fill(colorInventory, 0);
        
        /* Run GraphChi ... */
        GraphChiEngine<Integer[], Integer[]> engine = new GraphChiEngine<Integer[], Integer[]>(baseFilename, nShards);
        engine.setEdataConverter(new IntArrayConverter(numPartitions + 4));
        engine.setVertexDataConverter(new IntArrayConverter(2));
        engine.setModifiesInedges(false); // Important optimization
        engine.run(new JabeJaWeighted(), numberOfRounds);

        logger.info("Ready. Going to output...");

        /* Process output */
        PartitionAnalysis.writePartitions(logger, baseFilename, engine.numVertices(), engine.numVertices(), numPartitions);
        
        logger.info("Done!");
    }
}
