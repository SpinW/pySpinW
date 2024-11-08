# Architectural Decision Records

Significant architectural choices in the PySpinW project will use this workflow:

1. A pull request which creates a new markdown file in this folder with the proposed architectural design choice is opened.
2. Detailed design documentation supporting the decision will be placed separately into markdown files in the [design](../design) folder and should be referenced in the ADR file.
3. Project members and interested parties will comment on the pull request. If discussion meetings were convened about the decision, they should be minuted or summarised in comments on the PR.
4. The MD file is modified in light of the comments.
5. The PR is merged or closed without merging. Merging the PR signals that the proposed design decision is accepted.


## File name format

The files in this folder should named as `NNN-<short-title>` where `NNN` is a sequential number.


## Template

ADRs in PySpinW will use the template suggested in [this blog post](https://cognitect.com/blog/2011/11/15/documenting-architecture-decisions.html):

```
# Title

## Context

Detailed information should refer to design documents in the [design](../design) folder.

## Decision

"We will..."

## Status

Proposed / Accepted / Deprecated / Superseded

## Consequences
```

The ADR file should be short - not more than approximately one page if printed out. Detailed designs should be in the [design](../design) folder.
