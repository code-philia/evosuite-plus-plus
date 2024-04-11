package org.evosuite.coverage.branch;

import org.evosuite.TestGenerationContext;

public class BranchUtil {
    public static Branch getBranch(int branchId) {
        return BranchPool.getInstance(TestGenerationContext.getInstance().getClassLoaderForSUT())
                .getBranch(branchId);
    }

    public static int getLineOfBranchId(int branchId) {
        return getBranch(branchId).getInstruction().getLineNumber();
    }
}
