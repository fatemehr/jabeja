package edu.cmu.graphchi.apps.jabeja;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Iterator;
import java.util.logging.Logger;

import edu.cmu.graphchi.datablocks.IntArrayConverter;
import edu.cmu.graphchi.preprocessing.VertexIdTranslate;
import edu.cmu.graphchi.vertexdata.VertexAggregator;
import edu.cmu.graphchi.vertexdata.VertexIdValue;

public class PartitionAnalysis {

	/**
	 * Writes the computed partitions into separate output files 
	 * @author Fatemeh Rahimian
	 */
    
	public static void writePartitions(Logger logger, String baseFilename, int numShards, int numVertices, int numPartitions) {
		try {
	        Iterator<VertexIdValue<Integer[]>> iter = VertexAggregator.vertexIterator(numVertices, baseFilename,
	                new IntArrayConverter(2), VertexIdTranslate.identity());
	        
	        BufferedOutputStream []bos = new  BufferedOutputStream[numPartitions];
	        for (int i = 0; i < numPartitions; i++)
	        		bos[i] = new BufferedOutputStream(new FileOutputStream(baseFilename + "." + i));
	        
	        VertexIdValue<Integer[]> idVal;
	        while(iter.hasNext()) {
	            idVal = iter.next();
	            String s = idVal.getVertexId() + "\n";
	            bos[idVal.getValue()[1]].write(s.getBytes());
	        }
	        
	        for (int i = 0; i < numPartitions; i++) {
		        bos[i].flush();
		        bos[i].close();
	        }
	        
	        logger.info("PartitionAnalysis: output is written to files " + baseFilename + ".x, where x is an integer indicating the partition index");
		}
		catch (IOException e) {
			logger.info("PartitionAnalysis: can not write to file" + baseFilename + ".components\t" + e.getMessage());
		}
	}

}
