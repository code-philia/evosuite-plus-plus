package regression.objectconstruction.graphgeneration.testcase.multipath;

import java.io.File;
import java.lang.reflect.Method;
import java.util.Arrays;

import org.evosuite.Properties;
import org.evosuite.setup.DependencyAnalysis;
import org.evosuite.utils.MethodUtil;
import org.junit.Test;

import common.TestUtility;

public class MultiPathExampleTest {
	@Test
	public void testMultiPathExample1() {
		Class<?> clazz = regression.objectconstruction.graphgeneration.example.multipath.MultiPathExample.class;

		String methodName = "checkRules1";
		int parameterNum = 2;

		String targetClass = clazz.getCanonicalName();
		Method method = TestUtility.getTargetMethod(methodName, clazz, parameterNum);

		String targetMethod = method.getName() + MethodUtil.getSignature(method);
		
		Properties.TARGET_CLASS = targetClass;
		Properties.TARGET_METHOD = targetMethod;
		
		String cp = "target/test-classes";
		
		try {
			DependencyAnalysis.analyzeClass(Properties.TARGET_CLASS, Arrays.asList(cp.split(File.pathSeparator)));
		} catch (ClassNotFoundException | RuntimeException e) {
			e.printStackTrace();
		}
	}
	
	@Test
	public void testMultiPathExample2() {
		Class<?> clazz = regression.objectconstruction.graphgeneration.example.multipath.MultiPathExample.class;

		String methodName = "checkRules2";
		int parameterNum = 2;

		String targetClass = clazz.getCanonicalName();
		Method method = TestUtility.getTargetMethod(methodName, clazz, parameterNum);

		String targetMethod = method.getName() + MethodUtil.getSignature(method);
		
		Properties.TARGET_CLASS = targetClass;
		Properties.TARGET_METHOD = targetMethod;
		
		String cp = "target/test-classes";
		
		try {
			DependencyAnalysis.analyzeClass(Properties.TARGET_CLASS, Arrays.asList(cp.split(File.pathSeparator)));
		} catch (ClassNotFoundException | RuntimeException e) {
			e.printStackTrace();
		}
	}
}
