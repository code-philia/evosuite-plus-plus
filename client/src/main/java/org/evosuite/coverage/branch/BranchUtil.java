package org.evosuite.coverage.branch;

import org.evosuite.TestGenerationContext;
import org.evosuite.graphs.cfg.BytecodeInstruction;

public class BranchUtil {
    public static Branch getBranch(int branchId) {
        return BranchPool.getInstance(TestGenerationContext.getInstance().getClassLoaderForSUT())
                .getBranch(branchId);
    }

    public static int getLineOfBranchId(int branchId) {
        return getBranch(branchId).getInstruction().getLineNumber();
    }

    public static int getBranchBytecodeOffset(int branchId) {
        return getBranch(branchId).getInstruction().getBytecodeOffset();
    }

    public static int findBranchAtStatement(int statementId) {
        int countOfBranches = BranchPool.getInstance(TestGenerationContext.getInstance().getClassLoaderForSUT())
                .getAllBranchCount();
        for (int i = 0; i < countOfBranches; i++) {
            Branch branch = getBranch(i);
            if (branch.getInstruction().getBytecodeOffset() == statementId) {
                return i;
            }
        }
        return -1;
    }
}
