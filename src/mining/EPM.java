package mining;


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.BinomialDistribution;


public class EPM {
	/*
	 * Number of values of the class attribute
	 */
	private int m_classNum;
	/*
	 * Frequencies of each of the class labels
	 */
	private double[] m_classCounts;

	/*
	 * Minimum frequency threshold for each of the class labels
	 */
	private int[] m_minFrequencies;


	private double m_significanceThreshold; 
	/*
	 * Corresponding p-value for m_SignificanceThreshold, 
	 * which corresponds to the discriminative power threshold in the paper  
	 */
	private double m_pValueThreshold;

	private double m_minFrequencyThreshold;
	/*
	 * Maximum length of mined patterns
	 */
	private int m_maxSearchDepth;
	/* threshold used for information novelty
	 */
	private double m_binomial_p_threshold = 0.01;

	/*
	 * Separator used to separate columns in the input file
	 */
	public String m_saperator;



	/*
	 * The column index for the class attribute
	 */
	private int m_classAttIdx;
	/*


	/*
	 * Maximum number of class labels
	 */
	private int m_maxClassNum = 60;

	private boolean m_useBonferroniCorrection = true;
	/*
	 * Required for Bonferroni Correction
	 */
	private int[] m_testedPatternsCounts;


	private Map<Pair<Integer>, Double> m_likelihoodCache;

	/*
	 * Each feature-value (item) is stored as a unique index in FP-Tree
	 * m_indexToItem keeps the corresponding item of each of the indices  
	 */
	private Map<Integer, Item> m_indexToItem;

	private List<String> m_classList;


	public EPM(double minFrequency, int maxSearchDepth, 
			double significanceThreshold, int classIndex, String separator) {
		m_minFrequencyThreshold = minFrequency;
		m_classAttIdx = classIndex;
		m_maxSearchDepth = maxSearchDepth;
		m_saperator = separator;
		m_significanceThreshold = significanceThreshold;
		m_pValueThreshold = ChiSquare.chisqr(1, m_significanceThreshold);

		m_likelihoodCache = new HashMap<Pair<Integer>, Double>();

		m_classCounts = new double[m_maxClassNum];
		m_minFrequencies = new int[m_maxClassNum];

		m_indexToItem = new HashMap<Integer,Item>();
		m_classList = new ArrayList<String>();

		m_testedPatternsCounts = new int[m_maxSearchDepth];
	}


	/*
	 * The Item class represents a feature-value pair
	 */

	public class Item implements Serializable, Comparable<Item> {
		/** For serialization */
		private static final long serialVersionUID = -3372941834914147669L;

		/** The frequencies of the item in each of the class labels*/
		protected int[] m_counts;

		/** The name of the feature that the item corresponds to */
		protected String m_featureName;

		/** The corresponding index of the item*/
		protected int m_index;

		public Item (int index, String name, int[] counts) {
			m_index = index;
			m_featureName = name;
			m_counts = counts;
		}

		public Item (int index, String name) {
			m_index = index;
			m_featureName = name;
			m_counts = new int[m_classNum];
		}

		protected int getIndex(){
			return m_index;
		}


		private double getProbability(int targetClass){
			double sum = IntStream.of(m_counts).sum();
			return (m_counts[targetClass]/sum);
		}



		/**
		 * A string representation of this item.
		 * 
		 * @return a string representation of this item.
		 */
		public String toString() {
			return toString(false);
		}

		/**
		 * A string representation of this item.
		 * 
		 * @param printFrequency true if the frequency should be included.
		 * @return a string representation of this item. 
		 */
		public String toString(boolean printFrequency) {
			String freqStr = "";
			if (printFrequency){
				for (int i = 0; i < m_classNum; i++){
					freqStr += m_counts[i];
					if (i < m_classNum-1)
						freqStr += " , ";
				}
			}

			String result = m_featureName + (printFrequency? "("+freqStr+") " : "");

			return result;
		}

		public String shortDescription(){
			return m_featureName;
		}

		@Override
		public int compareTo(Item comp) {
			int thisSum =  IntStream.of(m_counts).sum();
			int compSum = IntStream.of(comp.m_counts).sum();

			if (thisSum == compSum) {
				return m_index - comp.m_index;
			}
			if (compSum < thisSum) {
				return -1;
			}
			return 1;
		}


		@Override
		public boolean equals(Object compareTo) {
			if (!(compareTo instanceof Item)) {
				return false;
			}

			Item b = (Item)compareTo;
			if (m_index == b.m_index) {
				for (int i = 0; i < m_classNum; i++)
					if (m_counts[i]!= b.m_counts[i])
						return false;
				return true;
			}

			return false;
		}

		@Override
		public int hashCode() {
			return (int) Math.pow(31 * m_index, IntStream.of(m_counts).sum()); 
		}
	}

	/**
	 * Class for maintaining a frequent item set, i.e. a pattern.
	 */

	public class Pattern 
	implements Serializable, Cloneable{

		private static final long serialVersionUID = 1L;

		/** List of items that are contained in this pattern */
		private List<Integer> m_itemIdxs;

		/** The frequency of this pattern in each of the class labels*/
		protected int[] m_patternClassCounts;

		/**The class label in which this pattern has the highest frequency */
		private int m_estimatedClass = -1;


		public Pattern() {
			m_itemIdxs = new ArrayList<Integer>();
			this.m_patternClassCounts = new int[m_classNum];
		}

		public Pattern(Pattern f) {
			m_itemIdxs = new ArrayList<Integer>();
			m_itemIdxs.addAll(f.m_itemIdxs);
			m_patternClassCounts = f.m_patternClassCounts;
			m_estimatedClass = f.m_estimatedClass;
		}


		public int getEstimatedClass(){
			if (m_estimatedClass == -1) {
				double max = m_patternClassCounts[0];
				m_estimatedClass = 0;
				for (int i = 1; i < m_classNum-1; i++){
					if (max < m_patternClassCounts[i]){
						max = m_patternClassCounts[i];
						m_estimatedClass = i;
					}
				}
			}

			return m_estimatedClass;
		}


		private double getProbability(int classIdx){
			double sum = IntStream.of(m_patternClassCounts).sum();
			return (m_patternClassCounts[classIdx]/sum);
		}

		private double getProbability(){
			double sum = IntStream.of(m_patternClassCounts).sum();
			return (m_patternClassCounts[getEstimatedClass()]/sum);
		}

		public void addItem(Integer idx) {
			m_itemIdxs.add(idx);
		}

		/**
		 * Set the frequency for this item set.
		 */
		public void setFrequency(int[] counts) {
			m_patternClassCounts = counts;
		}

		public int getFrequency(int classIdx) {
			return m_patternClassCounts[classIdx];
		}


		/**
		 * Computing the significance based on the contingency table of this pattern and the class label targetClass 
		 */
		private double getPatternSignificance(int targetClass){

			double patternCountInOtherClasses = 0;
			double otherClassesSum = 0;

			for (int i = 0; i < m_classNum; i++){
				if (i != targetClass){
					patternCountInOtherClasses+= m_patternClassCounts[i];
					otherClassesSum += m_classCounts[i];
				}
			}
			return logLikelihoodRatio(m_patternClassCounts[targetClass], patternCountInOtherClasses, m_classCounts[targetClass] - m_patternClassCounts[targetClass], otherClassesSum - patternCountInOtherClasses);
		}

		/**
		 * Similar to getPatternSignificance in which the targetClass is the class label in which this pattern has the highest frequency
		 */

		private double getPatternSignificance(){

			int targetClass = getEstimatedClass();
			double countInOtherClasses = 0;
			double otherClassesSum = 0;
			for (int i = 0; i < m_classNum; i++){
				if (i != targetClass){
					countInOtherClasses+= m_patternClassCounts[i];
					otherClassesSum += m_classCounts[i];
				}
			}
			return logLikelihoodRatio(m_patternClassCounts[targetClass], countInOtherClasses, m_classCounts[targetClass] - m_patternClassCounts[targetClass], otherClassesSum - countInOtherClasses);
		}


		public String toString(boolean displayStatistics) {
			String featureString =  " [";

			for (Integer itemIdx : m_itemIdxs){
				featureString += m_indexToItem.get(itemIdx).shortDescription() + " ";
			}
			featureString += "]("+ getEstimatedClass()+") ";

			if (displayStatistics){
				String countStr = "";
				for (int i = 0; i < m_classNum; i++)
					countStr += m_patternClassCounts[i] + " ";
				featureString += "(" + countStr +")";
			}
			return featureString;
		}

		@Override
		public String toString() {
			return toString(false);
		}

		public Object clone() {
			Pattern cloned = new Pattern(this);
			return cloned;
		}
	}

	class Pair<T> {
		T first;
		T second;

		public Pair(T first, T second){
			this.first = first;
			this.second = second;
		}

		@Override
		public boolean equals(Object o){
			if (o instanceof Pair) {
				@SuppressWarnings("unchecked")
				Pair<T> p = (Pair<T>) o;
				return first == p.first &&
						second == p.second;
			}
			return false;
		}

		@Override
		public String toString(){
			return first+"="+second;
		}
	}


	/**
	 * This class holds the counts for projected tree nodes
	 * and header lists.
	 */
	protected class ShadowCounts implements Serializable, Cloneable {

		/** For serialization */
		private final static long serialVersionUID = 4435433714185969155L;

		/** Holds the counts at different recursion levels */
		private Map<Integer, int[]> m_counts;

		public ShadowCounts(){
			m_counts = new HashMap<Integer, int[]>();
		}

		/**
		 * Get the count at the specified recursion depth.
		 * 
		 * @param recursionLevel the depth of the recursion.
		 * @return the count.
		 */
		public int[] getCount(int recursionLevel) {
			if (recursionLevel >= m_counts.size()) {
				return new int[m_classNum];
			} else {
				return m_counts.get(recursionLevel);
			}
		}


		/**
		 * Increase the count at a given recursion level.
		 * 
		 * @param recursionLevel the level at which to increase the count.
		 * @param incr the amount by which to increase the count.
		 */

		public void increaseCount(int recursionLevel, int[] incr){
			if (!m_counts.containsKey(recursionLevel)) {
				// new element
				int[] counts = new int[m_classNum];
				for (int i = 0; i < m_classNum; i++)
					counts[i] = incr[i];

				m_counts.put(recursionLevel, counts);
			} else {
				// otherwise increment the top
				int[] curr = m_counts.get(recursionLevel);

				for (int i = 0; i < m_classNum; i++)
					curr[i] += incr[i];
				m_counts.remove(recursionLevel);
				m_counts.put(recursionLevel, curr);
			}
		}

		/**
		 * Remove the count at the given recursion level.
		 * 
		 * @param recursionLevel the level at which to remove the count.
		 */
		public void removeCount(int recursionLevel) {
			m_counts.remove(recursionLevel);
		}

		@Override
		protected Object clone() throws CloneNotSupportedException {
			return super.clone();
		}
	}

	/**
	 * Calculate the Shannon entropy.
	 * @return The entropy value for the elements
	 */
	public double entropy(double... elements) {
		double sum = 0;
		for (double element : elements) {
			sum += element;
		}
		double result = 0.0;
		for (double x : elements) {
			if (x < 0) {
				throw new IllegalArgumentException("Should not have negative count for entropy computation: (" + x + ')');
			}
			int zeroFlag = (x == 0 ? 1 : 0);
			result += x * Math.log((x + zeroFlag) / sum);
		}
		return -result;
	}


	private double logLikelihoodRatio(double k11, double k12, 
			double k21, double k22) {
		Pair<Integer> p = new Pair<Integer>((int)k11, (int)k12);
		if (m_likelihoodCache.containsKey(p))
			return m_likelihoodCache.get(p);

		double rowEntropy = entropy(k11, k12) + entropy(k21, k22);
		double columnEntropy = entropy(k11, k21) + entropy(k12, k22);
		double matrixEntropy = entropy(k11, k12, k21, k22);
		if (rowEntropy + columnEntropy > matrixEntropy) {
			// round off error
			return 0.0;
		}
		double llr = 2.0 * (matrixEntropy - rowEntropy - columnEntropy);
		double formatted_llr = Double.parseDouble(String.format("%.2f", llr));
		m_likelihoodCache.put(p, formatted_llr);
		return formatted_llr;
	}




	/**
	 * FP-tree node.
	 * The implementation is mostly taken from that of Weka
	 */

	public class TreeNode implements Serializable {

		/** For serialization */
		private static final long serialVersionUID = 4396315323673737660L;

		/** link to the parent node */
		protected TreeNode m_parent;

		/** item at this node */
		protected int m_item;

		/** ID (for graphing the tree) */
		protected int m_ID;

		/** the children of this node */
		protected Map<Integer, TreeNode> m_children = 
				new HashMap<Integer, TreeNode>();

		/** counts associated with projected versions of this node */
		private ShadowCounts m_projectedInfo = new ShadowCounts();

		/**
		 * Construct a new node with the given parent link and item.
		 * 
		 * @param parent a pointer to the parent of this node.
		 * @param item the item at this node.
		 */
		public TreeNode(TreeNode parent, int item) {
			m_parent = parent;
			m_item = item;
		}

		/**
		 * Insert an item set into the tree at this node. Removes the first item
		 * from the supplied item set and makes a recursive call to insert the
		 * remaining items.
		 * 
		 * @param itemSet the item set to insert.
		 * @param headerTable the header table for the tree.
		 * @param incr the amount by which to increase counts.
		 */
		public void  addItemSet(Collection<Item> itemSet, 
				Map<Integer, Header> headerTable, int[] count) {
			Iterator<Item> i = itemSet.iterator();
			if (i.hasNext()) {
				Item first = i.next();
				TreeNode aChild;
				int index = first.getIndex();

				if (!m_children.containsKey(index)) {
					// not in the tree, so add it.
					aChild = new TreeNode(this, index);
					m_children.put(index, aChild);
					// update the header
					if (!headerTable.containsKey(index)) {
						headerTable.put(index, new Header());
					}


					// append new node to header list
					headerTable.get(index).addToList(aChild);
				} else {
					// get the appropriate child node
					aChild = m_children.get(index);
				}
				// update counts in header table
				headerTable.get(index).getProjectedCounts().increaseCount(0, count);
				// increase the child's count
				aChild.increaseProjectedCount(0, count);
				// proceed recursively
				itemSet.remove(first);
				aChild.addItemSet(itemSet, headerTable, count);
			}
		}

		/**
		 * Increase the projected count at the given recursion level at this
		 * node
		 * 
		 * @param recursionLevel the recursion level to increase the node count
		 * at.
		 * @param incr the amount by which to increase the count.
		 */
		public void increaseProjectedCount(int recursionLevel, int[] incr) {
			this.m_projectedInfo.increaseCount(recursionLevel, incr);
		}

		/**
		 * Remove the projected count at the given recursion level for this
		 * node.
		 * 
		 * @param recursionLevel the recursion level at which to remove the count.
		 */
		public void removeProjectedCount(int recursionLevel) {
			this.m_projectedInfo.removeCount(recursionLevel);
		}

		/**
		 * Get the projected count at the given recursion level for this node.
		 * 
		 * @param recursionLevel the recursion level at which to get the count.
		 * @return the count.
		 */
		public int[] getProjectedCount(int recursionLevel) {
			return m_projectedInfo.getCount(recursionLevel);
		}

		/**
		 * Get the parent node.
		 * 
		 * @return the parent node.
		 */
		public TreeNode getParent() {
			return m_parent;
		}

		/**
		 * Get the item at this node.
		 * 
		 * @return the item at this node.
		 */
		public int getItem() {
			return m_item;
		}    

		protected void extractCoOccuredAttributes(Set<Integer> partialPath, Map<Integer, Set<Integer>> res){

			if(!res.containsKey(m_item))
				res.put(m_item, new HashSet<Integer>());
			res.get(m_item).addAll(partialPath);

			if (m_children == null || m_children.size()==0)
				return;

			Set<Integer> newPath = new HashSet<Integer>();
			newPath.addAll(partialPath);
			newPath.add(m_item);

			for (TreeNode t : m_children.values()){
				t.extractCoOccuredAttributes(newPath, res);
			}
		}

		/**
		 * Return a textual description of this node for a given recursion
		 * level.
		 * 
		 * @param recursionLevel the recursion depth to use.
		 * @return a textual description of this node.
		 */
		public String toString(int recursionLevel) {
			return toString("", recursionLevel);
		}

		/**
		 * Return a textual description of this node for a given recursion
		 * level.
		 * 
		 * @param prefix a prefix string to prepend.
		 * @param recursionLevel the recursion level to use.
		 * @return a textual description of this node. 
		 */
		public String toString(String prefix, int recursionLevel) {
			StringBuilder buffer = new StringBuilder();
			buffer.append(prefix);
			buffer.append("|  ");
			buffer.append("((" +m_indexToItem.get(m_item).shortDescription());
			buffer.append(" (");
			buffer.append( getProjectedCount(0)[0] + " " + getProjectedCount(0)[1]);
			buffer.append(") \n");
			for (TreeNode node : m_children.values()) {
				buffer.append(node.toString(prefix + "|  ", recursionLevel));
			}
			return buffer.toString();
		}

		public String nodeString() {
			StringBuilder buffer = new StringBuilder();

			buffer.append("(" +m_indexToItem.get(m_item).toString());

			return buffer.toString();
		}

		public String nodeString(int recursionLevel) {
			StringBuilder buffer = new StringBuilder();
			buffer.append("(" +m_indexToItem.get(m_item).toString());

			return buffer.toString();
		}
	}

	boolean isFrequent(int[] count){
		for (int i = 0; i < m_classNum; i++)
			if (count[i] >= m_minFrequencies[i])
				return true;
		return false;
	}



	public class Header implements Serializable {
		/** For serialization */
		private static final long serialVersionUID = -6583156284891368909L;

		/** The list of pointers into the tree structure */
		protected List<TreeNode> m_headerList = new LinkedList<TreeNode>();

		/** Projected header counts for this entry */
		protected ShadowCounts m_projectedHeaderCounts = new ShadowCounts();

		/**
		 * Add a tree node into the list for this header entry.
		 * 
		 * @param toAdd the node to add.
		 */
		public void addToList(TreeNode toAdd) {
			m_headerList.add(toAdd);
		}

		/**
		 * Get the projected counts for this header entry.
		 * 
		 * @return the projected counts for this header entry.
		 */
		public ShadowCounts getProjectedCounts() {
			return m_projectedHeaderCounts;
		}
		public List<TreeNode> getHeaderList() {
			return m_headerList;
		}
	}

	/**
	 * Root of the FPTree
	 */

	public class TreeRoot extends TreeNode {

		/** For serialization */
		private final static long serialVersionUID = 632150939785333297L;

		/**
		 * Stores a header entry for an FPTree 
		 */

		/** Stores the header table as mapped Header entries */
		protected TreeMap<Integer, Header> m_headerTable = 
				new TreeMap<Integer, Header>();

		protected Map<Integer, Set<Integer>> m_coOccuredAtts =
				new HashMap<Integer, Set<Integer>>();

		/**
		 * Create a new FPTreeRoot.
		 */
		public TreeRoot() {
			super(null, -1);
		}

		/**
		 * Insert an item set into the tree.
		 * 
		 * @param itemSet the item set to insert into the tree.
		 * @param incr the increment by which to increase counters.
		 */
		public void addItemSet(Collection<Item> itemSet, int[] count) {
			super.addItemSet(itemSet, m_headerTable, count);
		}


		public TreeMap<Integer, Header> getHeaderTable() {
			return m_headerTable;
		}

		public boolean isEmpty(int recursionLevel, int class_index) {
			for (TreeNode c : m_children.values()) {
				for (int i = 0; i < m_classNum; i++)
					if (c.getProjectedCount(recursionLevel)[i] > 0) {
						return false;
					}
			}
			return true;
		}

		public String toString(String pad, int recursionLevel) {
			StringBuilder result = new StringBuilder();
			result.append(pad);
			result.append("+ ROOT\n");

			for (TreeNode node : m_children.values()) {
				result.append(node.toString(pad + "|  ", recursionLevel));
			}
			return result.toString();
		}

		/**
		 * Get a textual description of the header table for this tree.
		 * 
		 * @param recursionLevel the recursion level to use.
		 * @return a textual description of the header table for this
		 * tree at a given recursion level.
		 */
		public String printHeaderTable(int recursionLevel) {
			StringBuffer buffer = new StringBuffer();
			for (Integer itemIdx : m_headerTable.keySet()) {
				Item item = m_indexToItem.get(itemIdx);
				buffer.append(item.toString());
				buffer.append(" : ");
				int[] count = m_headerTable.get(itemIdx).getProjectedCounts().getCount(recursionLevel);
				String countStr = "";
				for (int i = 0; i < m_classNum; i++)
					countStr += count[i] + " " ;
				buffer.append( "("+ countStr+")");
				buffer.append("\n");
			}
			return buffer.toString();
		}


		protected void extractCoOccuredAttributes(){
			for (TreeNode t : m_children.values()){
				t.extractCoOccuredAttributes(new HashSet<Integer>(), m_coOccuredAtts);
			}
		}

		protected Set<Integer> getCoOccuredAtts(Integer att, Map<Integer, Header> headerTable, 
				int recursionLevel){
			Set<Integer> nextLevelAtts = new HashSet<Integer>();

			if (m_coOccuredAtts.containsKey(att))
				for (Integer index : m_coOccuredAtts.get(att)){
					int[] count = headerTable.get(index).getProjectedCounts().getCount(recursionLevel);
					if (isFrequent(count)){
						nextLevelAtts.add(index);
					}
				}
			return nextLevelAtts;
		}

	}

	/*
	 * Get frequent items of the data
	 */
	protected Map <String, Integer> getFrequentItems(String dataFile, double minFrequency) throws Exception {
		Map <String, Integer> singletons = new HashMap<String, Integer>();
		try {		
			BufferedReader reader = new BufferedReader(new FileReader(dataFile));

			String line;
			Map <String, int[]> itemCount = new HashMap<String, int[]>();

			while( ((line = reader.readLine())!= null)){ 

				if (line.isEmpty() == true || line.startsWith("@") || line.startsWith("#")) {
					continue;
				}
				// split the transaction into items
				String[] lineSplited = line.split(m_saperator);
				String classAtt = lineSplited[m_classAttIdx].trim();
				if (!m_classList.contains(classAtt)){
					m_classList.add(classAtt);
				}
				int classIdx = m_classList.indexOf(classAtt);
				m_classCounts[classIdx]++;

				for(int i = 0; i < lineSplited.length; i++){
					if (i != m_classAttIdx){
						String itemName = lineSplited[i].trim();
						if (!itemCount.containsKey(itemName))
							itemCount.put(itemName, new int[m_maxClassNum]);

						int[] prevCounts = itemCount.get(itemName);
						prevCounts[classIdx]++;
						itemCount.put(itemName, prevCounts);
					}
				}
			}

			m_classNum = m_classList.size();

			for (int i = 0; i < m_classNum; i++){
				m_minFrequencies[i] = (minFrequency >= 1)? (int)minFrequency : (int)Math.ceil(minFrequency * m_classCounts[i]);
			}

			int item_idx = 0;

			for (String item_name : itemCount.keySet()){
				if (isFrequent(itemCount.get(item_name))){
					Item item = new Item(item_idx, item_name, itemCount.get(item_name));
					m_indexToItem.put(item_idx, item);
					singletons.put(item_name, item_idx);
					item_idx++;
				}
			}

			reader.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
		return singletons;
	}

	public TreeRoot buildTree(Map <String, Integer> singletons, String dataFile) {
		TreeRoot tree = new TreeRoot();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(dataFile));

			String line;
			// for each line (transaction) until the end of file
			while( ((line = reader.readLine())!= null)){ 
				if (line.isEmpty() == true || line.startsWith("@") || line.startsWith("#")) {
					continue;
				}
				ArrayList<Item> transaction = new ArrayList<Item>();

				String[] lineSplited = line.split(m_saperator);
				String classAtt = lineSplited[m_classAttIdx].trim();
				int classIdx = m_classList.indexOf(classAtt);

				for(int i = 0; i < lineSplited.length; i++){
					if (i != m_classAttIdx){
						String item_name = lineSplited[i].trim();
						if (singletons.containsKey(item_name)) {
							transaction.add(m_indexToItem.get(singletons.get(item_name)));
						}
					}
				}

				if (transaction.size() > 0){
					Collections.sort(transaction);

					int[] count = new int[m_classNum];
					count[classIdx] = 1;

					tree.addItemSet(transaction, count);

					transaction = null;
				}
			}
			// close the input file
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return tree;
	}

	protected void extractInformativePatterns(TreeRoot tree, Integer baseItem, Pattern conditionalPattern,
			List<Pattern> informativePatterns, int recursionLevel) {
		TreeMap<Integer, Header> headerTable = tree.getHeaderTable();
		Header header = headerTable.get(baseItem);

		int[] frequency = header.getProjectedCounts().getCount(recursionLevel);

		Pattern newPattern = 
				(Pattern) conditionalPattern.clone();
		newPattern.addItem(baseItem);
		newPattern.setFrequency(frequency);

		if (isFrequent(frequency)){

			if (isInformative(newPattern, newPattern.m_itemIdxs.size())){
				informativePatterns.add(newPattern);
			}

			if (newPattern.m_itemIdxs.size() >= m_maxSearchDepth){
				return;
			}

			//Build the conditional tree
			for (TreeNode n : header.getHeaderList()) {
				// push count up path to root
				int[] currCount = n.getProjectedCount(recursionLevel);
				TreeNode temp = n.getParent();
				while (temp != tree) {
					// set/increase for the node
					temp.increaseProjectedCount(recursionLevel+1, currCount);

					// set/increase for the header table
					headerTable.get(temp.getItem()).getProjectedCounts().increaseCount(recursionLevel+1, currCount);

					temp = temp.getParent();
				}
			}

			Set<Integer> next_level_atts = tree.getCoOccuredAtts(baseItem, headerTable, recursionLevel+1);

			for (Integer next_candidate_att : next_level_atts){
				extractInformativePatterns(tree, next_candidate_att , newPattern, informativePatterns, recursionLevel+1);
			}
			// reverse the propagated counts
			for (TreeNode n : header.getHeaderList()) {
				TreeNode temp = n;
				while (temp != tree) {
					temp.removeProjectedCount(recursionLevel+1);							
					temp = temp.getParent();
				}
			}

			// reverse the propagated counts in the header list
			// at this recursion level
			for (Header h : headerTable.values()) {
				h.getProjectedCounts().removeCount(recursionLevel+1);
			} 

		}
	}

	private boolean isInformative(Pattern f, int level){
		int estimatedClass = f.getEstimatedClass();

		m_testedPatternsCounts[level-1]++;

		if (f.getPatternSignificance(estimatedClass) < m_significanceThreshold){
			return false;
		}

		int n = IntStream.of(f.m_patternClassCounts).sum();

		double maxP = -1;
		for (Integer item_idx : f.m_itemIdxs){
			double itemP = m_indexToItem.get(item_idx).getProbability(estimatedClass); 
			if (itemP > maxP)
				maxP = itemP;
		}

		if (maxP != -1 && level > 1){
			BinomialDistribution bd = new BinomialDistribution(n, maxP);
			int k = (int)f.m_patternClassCounts[estimatedClass];
			double p = (1-bd.cumulativeProbability(k));
			if (p >= m_binomial_p_threshold){
				return false;
			}
			else
				return true;
		}
		else
			return true;
	}


	public List<Pattern> mineAllFeatures(String dataFile) {

		List<Pattern> informativePatterns = new ArrayList<Pattern>();

		try {
			Map<String, Integer> frequentItems = getFrequentItems(dataFile, m_minFrequencyThreshold);
			System.out.println("Number of frequent items: " + frequentItems.size());
			TreeRoot tree = buildTree(frequentItems, dataFile);
			System.out.println("FP-Tree is built");

			System.out.println("Extracting cooccuring attribute-values");
			tree.extractCoOccuredAttributes();
			System.out.println("Start feature mining ...");

			/*
			 * For each frequent items I_j, mine the set of informative patterns that contain I_j 
			 */

			for (String item_name : frequentItems.keySet()){
				extractInformativePatterns(tree, frequentItems.get(item_name), new Pattern(), informativePatterns, 0);
			}


			return postPruning(informativePatterns);

		} catch (Exception e) {
			e.printStackTrace();
		}
		return informativePatterns;
	}

	private List<Pattern> postPruning(List<Pattern> inputFeatSet){

		double[] levelCriticalValues = new double[m_maxSearchDepth];

		if (m_useBonferroniCorrection){
			for (int i = 0; i < m_maxSearchDepth; i++){
				levelCriticalValues[i] = m_pValueThreshold / (Math.pow(2, i+1) * m_testedPatternsCounts[i]);
				if (i > 0)
					levelCriticalValues[i] = Math.min(levelCriticalValues[i], levelCriticalValues[i-1]);
			}
		}

		List<Pattern> prunedFeatures = new ArrayList<Pattern>();

		prunedFeatures.addAll(inputFeatSet);

		Set<Pattern> toBeRemoved = new HashSet<Pattern>();

		for (int i = 0; i < prunedFeatures.size()-1; i++){
			Pattern f = prunedFeatures.get(i);
			if (m_useBonferroniCorrection){
				int level = f.m_itemIdxs.size(); 
				if (ChiSquare.chisqr(1, f.getPatternSignificance()) > levelCriticalValues[level-1]){
					toBeRemoved.add(f);
					continue;
				}
			}

			for (int j = i+1; j < prunedFeatures.size(); j++){
				Pattern s = prunedFeatures.get(j);

				if (secondIsMoreSpecefic(f, s)){

					int sum = IntStream.of(s.m_patternClassCounts).sum();

					BinomialDistribution bd = new BinomialDistribution(sum, f.getProbability());
					int k = (int)s.m_patternClassCounts[s.getEstimatedClass()];
					double p = (1-bd.cumulativeProbability(k));
					if (p >= m_binomial_p_threshold){
						toBeRemoved.add(s);
					}						

				}
				else if (secondIsMoreSpecefic(s, f)){

					int sum = IntStream.of(f.m_patternClassCounts).sum(); 
					BinomialDistribution bd = new BinomialDistribution(sum, s.getProbability());
					int k = (int)f.m_patternClassCounts[f.getEstimatedClass()];
					double p = (1-bd.cumulativeProbability(k));
					if (p >= m_binomial_p_threshold){
						toBeRemoved.add(f);
					}
				}
			}
		}

		prunedFeatures.removeAll(toBeRemoved);

		return prunedFeatures;

	}

	private boolean secondIsMoreSpecefic(Pattern f1, Pattern f2){
		if (f2.m_itemIdxs.size() > f1.m_itemIdxs.size() && f2.m_itemIdxs.containsAll(f1.m_itemIdxs))
			return true;
		return false;

	}

	public Set<Set<String>> getStringRep (List<Pattern> features){
		Set<Set<String>> strRep = new HashSet<Set<String>>();

		for (Pattern f : features){
			Set<String> featureRep = new HashSet<String>();
			for (Integer item_idx : f.m_itemIdxs){
				Item item = m_indexToItem.get(item_idx);
				featureRep.add(item.m_featureName);
			}
			strRep.add(featureRep);
		}

		return strRep;
	}

	public static boolean isDouble(String s){
		String decimalPattern = "([0-9]*)\\.([0-9]*)";  
		return java.util.regex.Pattern.matches(decimalPattern, s);
	}

	public static boolean isInteger(String s){
		String integerPattern = "([0-9]+)";  
		return java.util.regex.Pattern.matches(integerPattern, s);
	}



	public static void main(String[] args){
		String trainFile, outputFile=null;
		int search_depth;
		double min_support;

		if (args.length < 3 || 
				(!(isDouble(args[1]) || isInteger(args[1])) && !isInteger(args[3]))){
			System.out.println("Invalid arguments. \n The correct usage is as follows:\n\n "
					+"java mining.EPM file-name minimum-frequency pattern-maximum-length [output-file]\n\n"
					+ " (1) file-name: the path to the feature-value file \n\n "
					+ "Each line of the file contains the list of feature-values for one sample and its corresponding class label (similar to libsvm data format).\n "
					+ "Each column contains one feature-value. The column separator is by default set to whitespace.\n "
					+ "Change the separator argument in the EPM function if other separator is used.\n "
					+ "By default, the first column is assigned for the class label.\n "
					+ "Change the classIndex argument in the EPM function, if the class label is included in another column.\n\n "

					+ "(2) minimum-frequency: the minimum frequency threshold\n\n "
					+ "If minimum-frequency>=1, the absolute value would be used. \n "
					+ "e.g. \"java mining.EPM file-path 2 3\" mines all patterns that occur more than two times in either of the classes.\n "
					+ "If minimum-frequency \\in (0,1), then the minimum frequency threshold is computed for each class separately;\n "
					+ "the minimum frequency threshold for class c_i would be (minimum-frequency*the number of samples that have c_i class label)\n\n "
					+ ""
					+ "(3) pattern-maximum-length: the threshold for the maximum length of patterns \n\n "
					+ "e.g. if pattern-maximum-length is set to 3, EPM will not mine any patterns with more than three feature-values\n\n "

					+ "(4) output-file (optional): the name of the file for saving the list of informative patterns"
					);

			System.exit(-1);
		}

		trainFile = args[0];
		min_support = Double.parseDouble(args[1]);
		search_depth = Integer.parseInt(args[2]);
		if (args.length > 3)
			outputFile = args[3];


		EPM extractor = new EPM(min_support, search_depth, 6.6, 0, "\\s+");
		long startTime = System.currentTimeMillis();
		List<Pattern> allFeatures = extractor.mineAllFeatures(trainFile);
		long time = (System.currentTimeMillis() - startTime)/1000;
		Set<Set<String>> stringRep = extractor.getStringRep(allFeatures);

		if (outputFile != null){

			try {
				PrintWriter pw = new PrintWriter(outputFile);

				for (Set<String> s: stringRep){
					pw.println(Arrays.toString(s.toArray()));
					pw.println("----------------");
				}
				
				pw.close();

			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}

		}

		System.out.println(stringRep.size() + " patterns mined in " + time + " second");

	}
}
