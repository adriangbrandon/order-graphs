def edges(start, s):
    out, stack, cur = [], [], start
    for c in s:
        if c.isdigit():
            last = c + cur.replace(c, "")
            out.append((cur, c, last))
        elif c == "[":
            stack.append(cur)
            cur = last
        elif c == "]":
            cur = stack.pop()
    return out
