Serialisation Methods
=====================

The data structures in pySpinW are not a simple tree, but have shared leaves/branches. That is to say, it **does not** look like this...

```mermaid
   graph TD;
      A==>B==>C==>D;
      B==>E;
      C==>F;
      C==>G;

```

But instead looks more like this...

```mermaid
   graph TD;
      A==>B==>C==>D;
      B==>E==>D;
      B==>F==>G;
```

A common strategy for serialisation for a tree like structure is to recursively serialise/deserialise the tree, starting at the root node.
First we call `serialise` on the root node, then the root node calls `serialise` on its children and collects it into an object, the children
then call `serialise` on their children, collecting the data, and so on.

```mermaid
   graph TD;
      A-->|serialise|B-->|serialise|C;
      A-->|serialise|D-->|serialise|E;
      C-->|json|B-->|json|A;
      E-->|json|D-->|json|A;

```


However, if we do this with the data structures in pySpinW, shared nodes will be duplicated. Instead we need to do something different. 
Whilst doing this we would like to keep to a tree-like serialisation much as possible, as we only deviate a small bit.

Contexts
========

The strategy used to serialise pySpinW data is to use "contexts". 
When we come to a node that can be a duplicate, instead of serilising diretly, we interact with a context object:
1) We serialise the data for our node, but don't return it
2) Store the data with an id in the context object
3) Instead of returning the serialised data, we return a reference containing the id to caller
4) At the end of the serialisation, the data in the context is serialised in parallel


```mermaid
   graph TD;
      A-->|serialise|B-->|serialise|D;
      A-->|serialise|C-->|serialise|D;
      D-->|json|context-->|reference|B-->|json|A;
      D-->|json|context-->|reference|C-->|json|A;
```

TODO: Better diagram
