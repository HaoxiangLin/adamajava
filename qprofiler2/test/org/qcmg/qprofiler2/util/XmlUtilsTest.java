package org.qcmg.qprofiler2.util;


import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.lang.reflect.ParameterizedType;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

//import static org.junit.jupiter.api.Assertions.assertAll;
import javax.xml.parsers.ParserConfigurationException;

import org.junit.Test;
import org.qcmg.common.util.QprofilerXmlUtils;
import org.w3c.dom.Element;

public class XmlUtilsTest {
		
	@Test
	public void testReadGroupElement() throws ParserConfigurationException {								
//		assertAll(  ()-> {
//			Element ele = XmlUtils.createReadGroupNode( QprofilerXmlUtils.createRootElement( XmlUtils.readGroupsEle, null) , "id" );
//			ele.getAttribute("RGID").equals("id");
//			ele.getNodeName().equals("readGroup");
//		} );
		
		try {
			
			Element ele = XmlUtils.createReadGroupNode( QprofilerXmlUtils.createRootElement( XmlUtils.readGroupsEle, null) , "id" );
			assertEquals(ele.getAttribute(XmlUtils.Sname), "id");
			assertEquals(ele.getNodeName(), "readGroup");

		}catch(ParserConfigurationException e) {
			fail("unexpected exception happened!");
		}
	}
	
	
	
	public void  createHeaderRecordsTest() {
		
		
		
	}
	
	
	
	@Test
	public void xuTest() {
		TreeMap<Integer, String> intStringMap = new TreeMap<>();		
		intStringMap.put(1, "one");
		
		
		maptest1(intStringMap );
		

	}
	
	
	private <T> void  maptest1(Map<T, String> tallys ) {
		
		T firstKey = tallys.keySet().stream().findFirst().get();
		if(firstKey instanceof Integer )
			System.out.println("it is Map<Integer,String>");
		
		if(firstKey instanceof Number )
			System.out.println("it is Map<Number,String>");	
		 
		if(firstKey instanceof String )
			System.out.println("it is Map<String,String>");	
		
	}
	
	
	private <T> void  maptest(Map<T, String> tallys, Class<T> type ) {
		if(type == String.class)
			System.out.println("string class");
		
		else if (type == Integer.class)
			System.out.println("integer class");
		
		
		
//		Class<T> clazz = (Class<T>) tallys.getClass() ;
//		System.out.println( "input is " + Arrays.toString(clazz.getTypeParameters()) );
//		System.out.println( "original is " + Arrays.toString(tallys.getClass().getTypeParameters()) );
//		System.out.println( "keyset class is  " +  tallys.keySet().getClass() );
//		
//		if(tallys.getClass() == new TreeMap<Integer,String>().getClass())
//			System.out.println("same");
		 
		
		
//	 	if(tallys.getClass() == new Map<Number,String>)
//	    	
//	    	if(clazz.getTypeParameters()) 
		
	}
	
	private <T> void  typetest( T  input ) {
		 
		System.out.println( "input is " +  input.getClass() );
 		 
		
		
//	 	if(tallys.getClass() == new Map<Number,String>)
//	    	
//	    	if(clazz.getTypeParameters()) 
		
	}
	
	
	public static Class<?> findSuperClassParameterType(Object instance, Class<?> classOfInterest, int parameterIndex) {
		  Class<?> subClass = instance.getClass();
		  while (subClass != subClass.getSuperclass()) {
		    // instance.getClass() is no subclass of classOfInterest or instance is a direct instance of classOfInterest
		    subClass = subClass.getSuperclass();
		    if (subClass == null) throw new IllegalArgumentException();
		  }
		  ParameterizedType parameterizedType = (ParameterizedType) subClass.getGenericSuperclass();
		  return (Class<?>) parameterizedType.getActualTypeArguments()[parameterIndex];
		}
}
