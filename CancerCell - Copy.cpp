#include "CancerCell.h"

CancerCell::CancerCell(double initialSize, double initialNutrientLevel, double initialMigrationSpeed,
    double nutrientThreshold, double maxCellSize, double baseGrowthRate,
    double timeStep, double divisionSizeMultiplier)
    : cellSize(initialSize), nutrientLevel(initialNutrientLevel), migrationSpeed(initialMigrationSpeed)
    , migrationDirection(0.0), nutrientThreshold(nutrientThreshold), maxCellSize(maxCellSize)
    , baseGrowthRate(baseGrowthRate), timeStep(timeStep), divisionSizeMultiplier(divisionSizeMultiplier) {
}



bool CancerCell::isNutrientAvailable() {

    if (nutrientLevel >= 0.0) {
        return true; 
    }
    else {
        return false;
    }
}

void CancerCell::increaseCellSize() {
    double nutrientFactor = nutrientLevel / nutrientThreshold;
    double sizeFactor = cellSize / maxCellSize;

    double growthRate = baseGrowthRate * nutrientFactor * (1 - sizeFactor);

    double deltaSize = growthRate * timeStep;
    cellSize += deltaSize;
    migrationSpeed *= cellSize / cellSize; 
}

double CancerCell::calculateNutrientConsumptionRate() {

    double maxConsumptionRate = 0.5; 
    double maxCellSize = 10.0;
    double consumptionRate = (cellSize / maxCellSize) * maxConsumptionRate;

    return consumptionRate;
}

void CancerCell::consumeNutrient() {

    double nutrientConsumptionRate = calculateNutrientConsumptionRate();
    double nutrientConsumed = nutrientConsumptionRate * cellSize;
    nutrientLevel -= nutrientConsumed;

    if (nutrientLevel < 0) {
        nutrientLevel = 0;
    }
}

bool CancerCell::shouldCellDivide() {

    bool condition1 = (cellSize >= cellSize * divisionSizeMultiplier);
    bool condition2 = (nutrientLevel >= nutrientThreshold);

    return (condition1 && condition2);
}


CancerCell* CancerCell::createNewCell() {

    CancerCell* newCell = new CancerCell(cellSize, nutrientLevel, migrationSpeed,
        nutrientThreshold, maxCellSize, baseGrowthRate,
        timeStep, divisionSizeMultiplier);

    return newCell;
}

void CancerCell::adjustCellProperties(CancerCell* newCell) {

    newCell->migrationSpeed = migrationSpeed;
}

void CancerCell::divideResources(CancerCell* newCell) {

    double dividedNutrient = nutrientLevel * 0.5;
    nutrientLevel -= dividedNutrient;
    newCell->nutrientLevel = dividedNutrient;
}

bool CancerCell::shouldCellMigrate() {


}

std::vector<Node*> CancerCell::getNeighboringNodes() {


}

void CancerCell::selectMigrationDirection(const std::vector<Node*>& neighboringNodes) {


}

void CancerCell::moveCellToNewPosition() {

}

void CancerCell::UpdateCellState() {
    if (isNutrientAvailable()) {
        increaseCellSize();
        consumeNutrient();
    }

    if (shouldCellDivide()) {
        CancerCell* newCell = createNewCell();
        adjustCellProperties(newCell);
        divideResources(newCell);
    }

    if (shouldCellMigrate()) {
        std::vector<Node*> neighboringNodes = getNeighboringNodes();
        selectMigrationDirection(neighboringNodes);
        moveCellToNewPosition();
    }

}
