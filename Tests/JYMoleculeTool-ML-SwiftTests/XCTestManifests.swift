import XCTest

#if !canImport(ObjectiveC)
public func allTests() -> [XCTestCaseEntry] {
    return [
        testCase(JYMoleculeTool_ML_SwiftTests.allTests),
    ]
}
#endif
