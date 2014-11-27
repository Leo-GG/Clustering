# ifndef LINK_H
# define LINK_H

/**
 * @file link.h
 * @brief Link class definition
 *
 * Defines the Link class and implements a few inline methods
 */

/**
 * @class Link
 * Represents a link between two Nodes
 * Holds two Nodes and the distance between them
 */
class Link
{
    private:

        float distance_;       // Distance between the two Nodes
        shared_ptr<Node> A_;    // Shared pointer to the first Node
        shared_ptr<Node> B_;    // Shared pointer to the second Node

    public:

        /**
        * Constructor
        * @param distance Distance between the two Nodes in the link
        * @param A Shared pointer to the first Node
        * @param B Shared pointer to the second Node
        */
        Link(shared_ptr<Node> const& A, shared_ptr<Node> const& B,
             float distance) : distance_(distance), A_(A), B_(B){};

        /**
        * Return the distance between the link Nodes
        * @return distance
        */
        float getDistance() { return distance_; };

        /**
        * Return a shared pointer to the first Node
        * @return A
        */
        shared_ptr<Node> getNodeA() { return A_; }

        /**
        * Return a shared pointer to the second Node
        * @return B
        */
        shared_ptr<Node> getNodeB() { return B_; }

};

# endif
