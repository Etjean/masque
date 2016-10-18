/*
 * RawHierarchyTree.java
 * 
 * Copyright 2006 Michigan State University Board of Trustees
 * 
 * Created on June 24, 2002, 2:36 PM
 */
package edu.msu.cme.rdp.classifier.train;

import edu.msu.cme.rdp.readseq.utils.orientation.GoodWordIterator;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;


/**
 * A RawHierarchyTree holds the raw taxon and sequences information.
 * A RawHierarchyTree can only have one type of children:
 *  child treenodes or child sequences.
 * @author  wangqion
 * @version
 */
public class RawHierarchyTree {

    private String name;
    private int leaveCount;		   // the number of sequence leavess directly belong to this treenode
    private RawHierarchyTree parent;
    private Map<String, RawHierarchyTree> subclasses = new HashMap();  // child treenodes
    private int[] wordOccurrence;  //size is 64k for word size 8
    private Taxonomy taxon;    // the unique id for this treenode
    private int genusIndex = -1;  // -1 means not a genus node
    //if this node is a genus node, the index of this genus in the genusNodeList
    private float avgCopyNumber = 0.0f; // the mean number of gene copies of the directly children of this node, default 0 means info not available

    /** Creates new RawHierarchyTree given the name and its parent. 
     * Note: a RawHierarchyTree can only have one type of children:
     *  child treenodes or child sequences.
     * Usually the lowest ranked nodes contain the sequence leaves.
     * */
    public RawHierarchyTree(String n, RawHierarchyTree p, Taxonomy tax) {
        name = n;
        taxon = tax;
        leaveCount = -1;
        addParent(p);
    }

    /** Adds the parent HierarchyTree, also adds this node to the parent tree as a child. */
    private void addParent(RawHierarchyTree p) {
        parent = p;
        if (parent != null) {
            parent.addSubclass(this);
        }
    }

    /** Adds a child treenode. */
    private void addSubclass(RawHierarchyTree c) {
        subclasses.put(c.getName(), c);
    }

    protected void setCopyNumber(float c){
        this.avgCopyNumber = c;
    }
    
    public float getCopyNumber(){
        return avgCopyNumber;
    }
    
    public boolean hasCopyNumber(){
        return (avgCopyNumber > 0);
    }
        
    /** Gets the name of the treenode. */
    public String getName() {
        return name;
    }

    /** Gets the parent treenode. */
    public RawHierarchyTree getParent() {
        return parent;
    }

    /** Gets the array of the subclasses if any. */
    public Collection<RawHierarchyTree> getSubclasses() {
        return subclasses.values();
    }

    /** Gets the child treenode with the given name. */
    public RawHierarchyTree getSubclassbyName(String n) {
        return subclasses.get(n);
    }

    /** Gets the size of the children,
     *  either taxon children or sequence leaves.
     */
    public int getSizeofChildren() {
        int size;
        if ((size = subclasses.size()) > 0) {
            return size;
        }
        return leaveCount;
    }

    /** Gets the size of the child treenodes. */
    public int getSizeofSubclasses() {
        return subclasses.size();
    }

    /** This method initiates the word occurrences from a sequence for the
     * lowest level of the hierarchy tree. Duplicate words from one sequence 
     * will be count only once.
     */
    public void initWordOccurrence(LineageSequence pSeq, float[] wordPriorArr) throws IOException {
        if (leaveCount < 0) {
            leaveCount = 1;
        } else {
            leaveCount++;
        }

        GoodWordIterator iterator = new GoodWordIterator(pSeq.getSeqString());
        if (wordOccurrence == null) {
            wordOccurrence = new int[iterator.getMask() + 1];
        }

        // create a temporary list and initialize the value to be -1;
        int[] wordList = new int[iterator.getNumofWords()];

        for (int i = 0; i < wordList.length; i++) {
            wordList[i] = -1;
        }

        int numUniqueWords = 0;  // indicate the number of unique words
        // duplicated words in one sequence are only counted once
        while (iterator.hasNext()) {
            int index = iterator.next();  // index is the actual integer representation of the word

            if (!isWordExist(wordList, index)) {
                wordList[numUniqueWords] = index;
                wordOccurrence[index]++;
                numUniqueWords++;
                // now add the word to the wordPriorArr
                // duplicated words in one sequence are only counted once
                wordPriorArr[index]++;
            }
        }
    }

    /** Checks if this word already been added to the wordOccurrence.
     * Returns true if found.
     */
    private boolean isWordExist(int[] wordList, int wordIndex) {
        for (int i = 0; i < wordList.length; i++) {
            if (wordList[i] == wordIndex) {
                return true;
            }
            if (wordList[i] == -1) {
                return false;
            }
        }
        return false;
    }

    /**
     * Returns the size of the array wordOccurrence.
     */
    public int getWordOccurrenceSize() {
        return wordOccurrence.length;
    }

    /** Gets the word occurrence for the given word index.
     */
    public int getWordOccurrence(int wordIndex) {
        return wordOccurrence[wordIndex];
    }

    /** Resets the array wordOccurrence to null.
     */
    public void releaseWordOccurrence() {
        wordOccurrence = null;
    }

    /** Counts the number of sequence leaves below this tree. */
    public int getLeaveCount() {
        if (!(leaveCount < 0)) {
            return leaveCount;
        }
        if (getSizeofSubclasses() <= 0) {
            return leaveCount;
        }
        leaveCount = 0;
        Iterator i = subclasses.values().iterator();
        while (i.hasNext()) {
            leaveCount += ((RawHierarchyTree) i.next()).getLeaveCount();
        }
        return leaveCount;
    }
    
    /* the number of lowest level nodes below */
    public int getGenusNodeCount(){
        if (getSizeofSubclasses() <= 0 || this.taxon.hierLevel.equalsIgnoreCase("GENUS")) {
            return 1;
        }
        int genusNodeCount = 0;
        Iterator i = subclasses.values().iterator();
        while (i.hasNext()) {
            genusNodeCount += ((RawHierarchyTree) i.next()).getGenusNodeCount();
        }
        return genusNodeCount;
    }
    
    public boolean isSingleton() {
        return (getLeaveCount() > 1) ? false : true;
    }

    /** Counts the number of non-singleton sequence leaves below this tree. */
    public int getNonSingletonLeaveCount() {

        if (isSingleton()) {
            return 0;
        }
        if (getSizeofSubclasses() <= 0) {
            return leaveCount;
        }
        int nonSingleton = 0;
        Iterator i = subclasses.values().iterator();
        while (i.hasNext()) {
            nonSingleton += ((RawHierarchyTree) i.next()).getNonSingletonLeaveCount();
        }
        return nonSingleton;
    }

    /** Returns the taxon object of this treenode.
     */
    public Taxonomy getTaxonomy() {
        return taxon;
    }

    /** Sets the genus index of this treenode.
     */
    public void setGenusIndex(int i) {
        genusIndex = i;
    }

    /** Returns the genus index of this treenode.
     */
    public int getGenusIndex() {
        return genusIndex;
    }
    
    /** get all the lowest level nodes in given hierarchy level starting from the given root
     */
    public void getNodeMap(String level, HashMap<String, RawHierarchyTree> nodeMap) {

        if (this.taxon.getHierLevel().equalsIgnoreCase(level)) {
            nodeMap.put(this.name, this);
            return;
        }
        //start from the root of the tree, get the subclasses.
        Collection al = new ArrayList();

        if ((al = this.getSubclasses()).isEmpty()) {
            return;
        }
        Iterator i = al.iterator();
        while (i.hasNext()) {
            ((RawHierarchyTree) i.next()).getNodeMap(level, nodeMap);
        }
    }
}
