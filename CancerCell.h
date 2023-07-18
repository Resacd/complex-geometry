#pragma once
#ifndef _CANCER_CELL_H_
#define _CANCER_CELL_H_

#include <vector>
#include "Node.h"

class CancerCell {
private:
    double cellSize;
    double nutrientLevel;
    double migrationSpeed;
    double migrationDirection;
    double nutrientThreshold;
    double maxCellSize;
    double baseGrowthRate;
    double timeStep;
    double divisionSizeMultiplier;
public:
    CancerCell(double initialSize, double initialNutrientLevel, double initialMigrationSpeed,
        double nutrientThreshold, double maxCellSize, double baseGrowthRate,
        double timeStep, double divisionSizeMultiplier);
    bool isNutrientAvailable();
    void increaseCellSize();
    double calculateNutrientConsumptionRate();
    void consumeNutrient();
    bool shouldCellDivide();
    CancerCell* createNewCell();
    void adjustCellProperties(CancerCell* newCell);
    void divideResources(CancerCell* newCell);
    bool shouldCellMigrate();
    std::vector<Node*> getNeighboringNodes();
    void selectMigrationDirection(const std::vector<Node*>& neighboringNodes);
    void moveCellToNewPosition();
    void UpdateCellState();
};

#endif // !_CANCER_CELL_H_
