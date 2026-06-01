#version 450 core

layout(triangles) in;
layout(triangle_strip, max_vertices = 12) out;

uniform vec2 viewportSize;
uniform float lineWidth;

out vec2 edgeUV;


void emitEdgeQuad(vec4 clipA, vec4 clipB)
{
    vec2 ndcA = clipA.xy / clipA.w;
    vec2 ndcB = clipB.xy / clipB.w;

    vec2 screenA = (ndcA * 0.5 + 0.5) * viewportSize;
    vec2 screenB = (ndcB * 0.5 + 0.5) * viewportSize;

    vec2 dir = normalize(screenB - screenA);
    vec2 perp = vec2(-dir.y, dir.x);

    vec2 offset = perp * (lineWidth * 0.5);

    vec2 p0 = screenA - offset;
    vec2 p1 = screenA + offset;
    vec2 p2 = screenB - offset;
    vec2 p3 = screenB + offset;

    vec4 v[4];

    v[0] = vec4(
        ((p0 / viewportSize) * 2.0 - 1.0) * clipA.w,
        clipA.z,
        clipA.w);

    v[1] = vec4(
        ((p1 / viewportSize) * 2.0 - 1.0) * clipA.w,
        clipA.z,
        clipA.w);

    v[2] = vec4(
        ((p2 / viewportSize) * 2.0 - 1.0) * clipB.w,
        clipB.z,
        clipB.w);

    v[3] = vec4(
        ((p3 / viewportSize) * 2.0 - 1.0) * clipB.w,
        clipB.z,
        clipB.w);

    gl_Position = v[0];
    edgeUV = vec2(0,0);
    EmitVertex();

    gl_Position = v[1];
    edgeUV = vec2(0,1);
    EmitVertex();

    gl_Position = v[2];
    edgeUV = vec2(1,0);
    EmitVertex();

    gl_Position = v[3];
    edgeUV = vec2(1,1);
    EmitVertex();

    EndPrimitive();
}

void main()
{
    emitEdgeQuad(gl_in[0].gl_Position,
                 gl_in[1].gl_Position);

    emitEdgeQuad(gl_in[1].gl_Position,
                 gl_in[2].gl_Position);

    emitEdgeQuad(gl_in[2].gl_Position,
                 gl_in[0].gl_Position);
}