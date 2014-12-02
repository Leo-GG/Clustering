/**
 * @file link_comparator.h
 * @brief LinkComparator class definition
 *
 * Defines the LinkComparator class and implements a few inline methods
 */


# ifndef LINK_COMPARATOR_H
# define LINK_COMPARATOR_H

/**
 * @class LinkComparator
 * A simple comparator for Links between Nodes
 */
class LinkComparator {

public:

    /**
     * Comparator method. Returns TRUE if the distance in the first link is
     * larger or equal to the distance in the second link
     * @param linkA first link to compare
     * @param linkB second link to compare
     * @return TRUE if distance in linkA>=distance in link B
     */
    bool operator()(Link& linkA, Link& linkB)
    {
        bool compare = (linkA.getDistance() >= linkB.getDistance()) ? true : false ;
        return compare;
    }
};

#endif
