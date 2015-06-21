package org.evosuite.strategy;

import java.util.ArrayList;
import java.util.List;

import org.evosuite.ProgressMonitor;
import org.evosuite.Properties;
import org.evosuite.coverage.FitnessFunctions;
import org.evosuite.coverage.TestFitnessFactory;
import org.evosuite.rmi.ClientServices;
import org.evosuite.ga.stoppingconditions.GlobalTimeStoppingCondition;
import org.evosuite.ga.stoppingconditions.MaxFitnessEvaluationsStoppingCondition;
import org.evosuite.ga.stoppingconditions.MaxGenerationStoppingCondition;
import org.evosuite.ga.stoppingconditions.MaxStatementsStoppingCondition;
import org.evosuite.ga.stoppingconditions.MaxTestsStoppingCondition;
import org.evosuite.ga.stoppingconditions.MaxTimeStoppingCondition;
import org.evosuite.ga.stoppingconditions.StoppingCondition;
import org.evosuite.ga.stoppingconditions.ZeroFitnessStoppingCondition;
import org.evosuite.statistics.RuntimeVariable;
import org.evosuite.testcase.TestFitnessFunction;
import org.evosuite.testsuite.TestSuiteChromosome;
import org.evosuite.testsuite.TestSuiteFitnessFunction;

/**
 * This is the abstract superclass of all techniques to generate a set of tests
 * for a target class, which does not neccessarily require the use of a GA.
 * 
 * Postprocessing is not done as part of the test generation strategy.
 * 
 * @author gordon
 *
 */
public abstract class TestGenerationStrategy {

	/**
	 * Generate a set of tests; assume that all analyses are already completed
	 * @return
	 */
	public abstract TestSuiteChromosome generateTests();
	
	/** There should only be one */
	protected final ProgressMonitor progressMonitor = new ProgressMonitor();

	/** There should only be one */
	protected ZeroFitnessStoppingCondition zeroFitness = new ZeroFitnessStoppingCondition();
	
	/** There should only be one */
	protected StoppingCondition globalTime = new GlobalTimeStoppingCondition();

    protected void sendExecutionStatistics() {
        ClientServices.getInstance().getClientNode().trackOutputVariable(RuntimeVariable.Statements_Executed, MaxStatementsStoppingCondition.getNumExecutedStatements());
        ClientServices.getInstance().getClientNode().trackOutputVariable(RuntimeVariable.Tests_Executed, MaxTestsStoppingCondition.getNumExecutedTests());
    }
    
    /**
     * Convert criterion names to test suite fitness functions
     * @return
     */
	protected List<TestSuiteFitnessFunction> getFitnessFunctions() {
	    List<TestSuiteFitnessFunction> ffs = new ArrayList<TestSuiteFitnessFunction>();
	    for (int i = 0; i < Properties.CRITERION.length; i++) {
	        ffs.add(FitnessFunctions.getFitnessFunction(Properties.CRITERION[i]));
	    }

		return ffs;
	}
	
	/**
	 * Convert criterion names to factories for test case fitness functions
	 * @return
	 */
	public static List<TestFitnessFactory<? extends TestFitnessFunction>> getFitnessFactories() {
	    List<TestFitnessFactory<? extends TestFitnessFunction>> goalsFactory = new ArrayList<TestFitnessFactory<? extends TestFitnessFunction>>();
	    for (int i = 0; i < Properties.CRITERION.length; i++) {
	        goalsFactory.add(FitnessFunctions.getFitnessFactory(Properties.CRITERION[i]));
	    }

		return goalsFactory;
	}
	
	/**
	 * Check if the budget has been used up. The GA will do this check
	 * on its own, but other strategies (e.g. random) may depend on this function.
	 * 
	 * @param chromosome
	 * @param stoppingCondition
	 * @return
	 */
	protected boolean isFinished(TestSuiteChromosome chromosome, StoppingCondition stoppingCondition) {
		if (stoppingCondition.isFinished())
			return true;

		if (Properties.STOP_ZERO) {
			if (chromosome.getFitness() == 0.0)
				return true;
		}

		if (!(stoppingCondition instanceof MaxTimeStoppingCondition)) {
			if (globalTime.isFinished())
				return true;
		}

		return false;
	}
	
	/**
	 * Convert property to actual stopping condition
	 * @return
	 */
	protected StoppingCondition getStoppingCondition() {
		switch (Properties.STOPPING_CONDITION) {
		case MAXGENERATIONS:
			return new MaxGenerationStoppingCondition();
		case MAXFITNESSEVALUATIONS:
			return new MaxFitnessEvaluationsStoppingCondition();
		case MAXTIME:
			return new MaxTimeStoppingCondition();
		case MAXTESTS:
			return new MaxTestsStoppingCondition();
		case MAXSTATEMENTS:
			return new MaxStatementsStoppingCondition();
		default:
			return new MaxGenerationStoppingCondition();
		}
	}

}