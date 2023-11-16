#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#define kase(tipo,stmt) case(tipo):{stmt;break;}

char *a[10] = {
    "in particular",
    "on the other hand",
    "however",
    "similarly",
    "in this regard",
    "as a resultant implication",
    "based on integral subsystem considerations",
    "for example",
    "thus",
    "in respect to specific goals"},

    *b[10] = {
    "a large portion of the interface coordinated communication",
    "a constant flow of effective information",
    "the characterization of specific criteria",
    "initiation of critical subsystem development",
    "the fully integrated test program",
    "the product configuration baseline",
    "any associated supporting element",
    "the incorporation of additional mission constraints",
    "the independent functional principle",
    "a primary interrelationship between system and/or subsystem technologies"},

    *c[10] = {
    "must utilize and be functionally interwoven with",
    "maximizes the probability of project success and minimizes the cost and time required for",
    "adds explicit performance limits to",
    "necessitates that urgent consideration be applied to",
    "requires considerable systems analysis and trade off studies to arrive at",
    "is further compounded when taking into account",
    "presents extremely interesting challenges to",
    "recognizes the importance of other systems and the necessity for",
    "effects a significant implementation of",
    "adds overriding performance constraints to"},

    *d[10] = {
    "the sophisticated hardware",
    "the anticipated next generation equipment",
    "the subsystem compatibility testing",
    "the structural design based on system engineering concepts",
    "the preliminary qualification limits",
    "the evolution of specification over a given time period",
    "the philosophy of commonality and standardization",
    "the top-down development method",
    "any discrete configuration mode",
    "the total system rationale"}; /* orders: abcd, dacb, bacd, adcb */

main()
{
    int n, order, w, x, y, z;

    srand(time(NULL));
    for (n = 0; n < 1000; n++)
    {
        if (!(n % 10)) printf("\n");
        w = rand() % 10;
        x = rand() % 10;
        y = rand() % 10;
        z = rand() % 10;
        order = rand() % 4;
        switch (order)
        {
            case 0:
                printf(" %c%s, %s %s %s.", a[w][0] & 0xDF, a[w] + 1, b[x], c[y], d[z]);
                break;
            case 1:
                printf(" %c%s, %s, %s %s.", d[w][0] & 0xDF, d[w] + 1, a[x], c[y], b[z]);
                break;
            case 2:
                printf(" %c%s, %s, %s %s.", b[w][0] & 0xDF, b[w] + 1, a[x], c[y], d[z]);
                break;
            case 3:
                printf(" %c%s, %s %s %s.", a[w][0] & 0xDF, a[w] + 1, d[x], c[y], b[z]);
                break;
        }
    }
}
