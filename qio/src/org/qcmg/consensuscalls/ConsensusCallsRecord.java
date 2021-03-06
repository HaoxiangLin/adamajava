/**
 * © Copyright The University of Queensland 2010-2014.  This code is released under the terms outlined in the included LICENSE file.
 */
//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, vhudson-jaxb-ri-2.2-147 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2013.10.25 at 10:52:20 AM EST 
//


package org.qcmg.consensuscalls;

import java.util.ArrayList;
import java.util.List;
public class ConsensusCallsRecord {

    protected String chr;
    protected int position;
    protected String alleleDiColor1;
    protected String alleleDiColor2;
    protected String reference;
    protected String genotype;
    protected double pValue;
    protected List<ConsensusCallsFlag> flag;
    protected int coverage;
    protected int nCounts1StAllele;
    protected int nCountsReferenceAllele;
    protected int nCountsNonReferenceAllele;
    protected int refAvgQV;
    protected int novelAvgQV;
    protected int heterozygous;
    protected String algorithm;
    protected String algorithmName;

    /**
     * Gets the value of the chr property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getChr() {
        return chr;
    }

    /**
     * Sets the value of the chr property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setChr(String value) {
        this.chr = value;
    }

    /**
     * Gets the value of the position property.
     * 
     */
    public int getPosition() {
        return position;
    }

    /**
     * Sets the value of the position property.
     * 
     */
    public void setPosition(int value) {
        this.position = value;
    }

    /**
     * Gets the value of the alleleDiColor1 property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getAlleleDiColor1() {
        return alleleDiColor1;
    }

    /**
     * Sets the value of the alleleDiColor1 property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setAlleleDiColor1(String value) {
        this.alleleDiColor1 = value;
    }

    /**
     * Gets the value of the alleleDiColor2 property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getAlleleDiColor2() {
        return alleleDiColor2;
    }

    /**
     * Sets the value of the alleleDiColor2 property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setAlleleDiColor2(String value) {
        this.alleleDiColor2 = value;
    }

    /**
     * Gets the value of the reference property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getReference() {
        return reference;
    }

    /**
     * Sets the value of the reference property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setReference(String value) {
        this.reference = value;
    }

    /**
     * Gets the value of the genotype property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getGenotype() {
        return genotype;
    }

    /**
     * Sets the value of the genotype property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setGenotype(String value) {
        this.genotype = value;
    }

    /**
     * Gets the value of the pValue property.
     * 
     */
    public double getPValue() {
        return pValue;
    }

    /**
     * Sets the value of the pValue property.
     * 
     */
    public void setPValue(double value) {
        this.pValue = value;
    }

    /**
     * Gets the value of the flag property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the flag property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getFlag().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link ConsensusCallsFlag }
     * 
     * 
     */
    public List<ConsensusCallsFlag> getFlag() {
        if (flag == null) {
            flag = new ArrayList<ConsensusCallsFlag>();
        }
        return this.flag;
    }

    /**
     * Gets the value of the coverage property.
     * 
     */
    public int getCoverage() {
        return coverage;
    }

    /**
     * Sets the value of the coverage property.
     * 
     */
    public void setCoverage(int value) {
        this.coverage = value;
    }

    /**
     * Gets the value of the nCounts1StAllele property.
     * 
     */
    public int getNCounts1StAllele() {
        return nCounts1StAllele;
    }

    /**
     * Sets the value of the nCounts1StAllele property.
     * 
     */
    public void setNCounts1StAllele(int value) {
        this.nCounts1StAllele = value;
    }

    /**
     * Gets the value of the nCountsReferenceAllele property.
     * 
     */
    public int getNCountsReferenceAllele() {
        return nCountsReferenceAllele;
    }

    /**
     * Sets the value of the nCountsReferenceAllele property.
     * 
     */
    public void setNCountsReferenceAllele(int value) {
        this.nCountsReferenceAllele = value;
    }

    /**
     * Gets the value of the nCountsNonReferenceAllele property.
     * 
     */
    public int getNCountsNonReferenceAllele() {
        return nCountsNonReferenceAllele;
    }

    /**
     * Sets the value of the nCountsNonReferenceAllele property.
     * 
     */
    public void setNCountsNonReferenceAllele(int value) {
        this.nCountsNonReferenceAllele = value;
    }

    /**
     * Gets the value of the refAvgQV property.
     * 
     */
    public int getRefAvgQV() {
        return refAvgQV;
    }

    /**
     * Sets the value of the refAvgQV property.
     * 
     */
    public void setRefAvgQV(int value) {
        this.refAvgQV = value;
    }

    /**
     * Gets the value of the novelAvgQV property.
     * 
     */
    public int getNovelAvgQV() {
        return novelAvgQV;
    }

    /**
     * Sets the value of the novelAvgQV property.
     * 
     */
    public void setNovelAvgQV(int value) {
        this.novelAvgQV = value;
    }

    /**
     * Gets the value of the heterozygous property.
     * 
     */
    public int getHeterozygous() {
        return heterozygous;
    }

    /**
     * Sets the value of the heterozygous property.
     * 
     */
    public void setHeterozygous(int value) {
        this.heterozygous = value;
    }

    /**
     * Gets the value of the algorithm property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getAlgorithm() {
        return algorithm;
    }

    /**
     * Sets the value of the algorithm property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setAlgorithm(String value) {
        this.algorithm = value;
    }

    /**
     * Gets the value of the algorithmName property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getAlgorithmName() {
        return algorithmName;
    }

    /**
     * Sets the value of the algorithmName property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setAlgorithmName(String value) {
        this.algorithmName = value;
    }

}
