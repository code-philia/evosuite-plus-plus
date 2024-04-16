package org.evosuite.coverage.branch;

import org.evosuite.testcase.TestChromosome;
import org.evosuite.testsuite.TestSuiteChromosome;

public class BranchCoverageCase implements java.io.Serializable {
    public Branch branch;
    public int branchId;
    public String className;
    public String methodName;
    public int line;
    public TestSuiteChromosome bestCoveredTestSuite;
    public int bestCoveredTestSuiteId;
    public TestChromosome bestCoveredTestCase;
    public int bestCoveredTestCaseId;
    public double distance;

    public BranchCoverageCase(Branch branch) {
        this.branch = branch;
        this.branchId = branch.getActualBranchId();
        this.className = branch.getClassName();
        this.methodName = branch.getMethodName();
        this.line = branch.getInstruction().getLineNumber();
        this.bestCoveredTestSuite = null;
        this.bestCoveredTestSuiteId = -1;
        this.bestCoveredTestCase = null;
        this.bestCoveredTestCaseId = -1;
        this.distance = Double.MAX_VALUE;
        
    }

    public void setBest(TestSuiteChromosome testSuite, int testSuiteId, TestChromosome testCase, int testCaseId,
            double distance) {
        this.bestCoveredTestSuite = testSuite;
        this.bestCoveredTestSuiteId = testSuiteId;
        this.bestCoveredTestCase = testCase;
        this.bestCoveredTestCaseId = testCaseId;
        this.distance = distance;
    }
    
    public String toString() {
        return "BranchCoverageCase [branch=" + branch + ", branchId=" + branchId + ", className=" + className
                + ", methodName=" + methodName + ", line=" + line
                + ", bestCoveredTestSuiteId=" + bestCoveredTestSuiteId
                + ", bestCoveredTestCaseId=" + bestCoveredTestCaseId + ", distance=" + distance + "]";
    }
}
